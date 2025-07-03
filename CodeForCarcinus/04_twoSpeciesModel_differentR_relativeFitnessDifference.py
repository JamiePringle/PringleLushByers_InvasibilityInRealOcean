from pylab import *
from numpy import *
import zarr
import xarray as xr
from collections import Counter
import time
import os
import multiprocessing as mp
from multiprocessing import shared_memory
import createLinearModel_module as cLM
from numba import jit,njit
import tqdm

# This code models the stochastic population dynamics of two species
# which may have different larval production R. This is used to model
# a difference in relative fitness. R is the number of larvae produced
# -- some of these larvae may be lost because they do not return to
# suitable habitat.
#
# The question that this code intends to answer is what is the
# likelihood that an introduction into a single region will persist
# over some time. This experiment is repeated for all the regions that
# we define in '01_makeInitalIntroductionRanges.py', so we can
# understand how this likelihood differes between regions. Unlike the
# neutral case, it is necessary to run the model seperately for each
# introduction in each region.
#
# For speed, the conectivity matrices are converted from ones where
# location is defined by the (nx,ny) grid locations in the underlying
# circulation model to a linear domain defined by a single index. This
# is done for speed, and does not affect anything -- nothing in the
# underlying connectivity of the habitats is changed.

#assumes the model output directory modelOutputNeutral exists
assert os.path.exists('modelOutputRelativeFitness'), 'please create directory modelOutputRelativeFitness to store model output'  


#initialize random number generator
rng=np.random.default_rng()

#we load the large connectivity data in global scope, so we don't have to pass it into the
#parallel functions each time we call the. 
#print('    Loading model',flush=True)
#load a connectivity matrix and convert to linear models
ConnectivityModelName='E_CmaenasHab_depth1_minPLD40_maxPLD40_months5_to_6'
EfileName='transposes/'+ConnectivityModelName+'.zarr'
EwhereTo,EnumTo,EfracReturn,Ecumsum,\
    nxny2nlin,nlin2nxny,nxny2lonLat,nlin2lonLat=cLM.makeLinearModel(EfileName)
#print('       done loading model')

#this is a function which quickly counts how many of each value in
#whereSettleInDomain there are. It uses Numba to run faster, and
#is faster than the Counter class. 
@njit
def numbaFastCount(whereSettleInDomain,thisWhereSettle):
    for n in whereSettleInDomain:
        thisWhereSettle[n]+=1
    return None

#for better profiling, and future parallelization, lets move the computational core to a function
#@profile
def findWhereSettle(P,lowerBound,upperBound,R0,R1,Ndomain,Nspecies):
    #see below for what these parameters are. 

    tic=time.time()

    whereSettle=zeros((Ndomain,Nspecies),dtype=int)

    #only work on species that exist
    Psum=P.sum(axis=0)
    
    #now loop over species
    for nsp in range(lowerBound,upperBound):
        if Psum[nsp]>0:
            if nsp==0:
                R=R0
            elif nsp==1:
                R=R1
            else:
                assert False,'only set up for two species'

            #What is total number of larvae launched from each point which survive
            totalSurvive=R*P[:,nsp]*EfracReturn #note, can be much less than 1!

            #now convert totalSurvive into an array of integers. The
            #likelyhood that it is rounded up is the fractional part, the
            #likelyhood that it is round down is the 1-fractional part. We
            #do this by adding to totalSurvive a random uniform number
            #between 0 and 1, and then rounding down.
            intSurvive=floor(rng.random(Ndomain)+totalSurvive).astype(int)

            thisWhereSettle=zeros((Ndomain,),dtype=int)
            for n in arange(Ndomain)[P[:,nsp]>0]:
                #make array of integers which is list of which entries in EwhereTo of where the larvae are going
                #for future optimization, the rng.random() takes about 1/3 of the time of the whereGoList.
                #this is the main hotspot of the code, and takes about 60% of the total time
                whereGoList=searchsorted(Ecumsum[n],rng.random(intSurvive[n]))

                #now make array of what indices in the linear domain the larvae are going into
                whereSettleInDomain=EwhereTo[n][whereGoList]

                #where settle contains the number of larve from the species which reach
                #each point in the domain.
                if False:
                    #this is actually much slower, by a factor of two... But is this true
                    #for all model configurations? Perhaps experiment and see for your case
                    thisWhereSettle=thisWhereSettle+bincount(whereSettleInDomain,minlength=Ndomain)
                elif False:
                    for n in whereSettleInDomain:
                        thisWhereSettle[n]+=1
                elif False:
                    jnk=Counter(whereSettleInDomain)
                    for j in jnk:
                        thisWhereSettle[j]+=jnk[j]
                else:
                    numbaFastCount(whereSettleInDomain,thisWhereSettle)
                #breakpoint()

            #put answer into whereSettle
            whereSettle[:,nsp]=thisWhereSettle

    #print('      done with a processes in mp.Pool in',time.time()-tic,flush=True)
    return whereSettle

#write a function that determines how many of the larvae in
#whereSettle survive to be adults in P.  now there are two
#possibilities. Where the total number of larvae reaching a location
#is less than or equal to Pmax, they all survive. Where it is greater,
#choose Pmax survivors randomly from the larvae. Should make sure
#number of recruits in each species is not greater than number of
#larvae which reach (really? or can the average number be right?). Do
#this in a function so it is easy to profile. NOTE WELL, assume all
#fitness differences are in the fecundity and dispersal, and thus
#handled by findWhereSettle()
#@profile
def whichLarvaeSurvive(whereSettle,lowerBound,upperBound,Pmax,Ndomain,Nspecies):
    #inputs:
    #  P: population matrix (Ndomain,Nspecies). Pass in name of shared memory in which it is stored
    #  whereSettle: (Ndomain,Nspecies) number of larvae of each species
    #     which reach each locations. Pass in name of shared memory in which it is stored
    #  Pmax, Ndomain, and Nspecies as defined below

    P=zeros((Ndomain,Nspecies),dtype=int)

    #find total number of recruits, and put them in population
    #if the number is less or equal than Pmax
    totalRecruits=whereSettle.sum(axis=1)

    #find which regions this call to whichLarvaeSurvive() will deal with...
    #make an indx of all False's, then fill area we are interested in with True
    indxWorryAbout=zeros((Ndomain,),dtype=bool)
    indxWorryAbout[lowerBound:upperBound]=True
    
    #don't worry about empty habitat? Does this save time? Then have to make a new index below...
    indx=logical_and(totalRecruits<=Pmax,totalRecruits>0)
    indx=logical_and(indx,indxWorryAbout) #only process between lower and upperBound
    P[indx,:]=whereSettle[indx,:]

    #now the hard part. if totalRecruits>Pmax, then choose randomly
    #which recruits survive. How to do this quickly?
    indx=totalRecruits>Pmax
    indx=logical_and(indx,indxWorryAbout) #only process between lower and upperBound
    whichPlaces=arange(Ndomain,dtype=int)[indx]
    for nin in whichPlaces:
        #cumulative sum of fraction of settlers by species
        thisCumsum=cumsum(whereSettle[nin,:])/totalRecruits[nin]

        #randomly choose Pmax survivors by fraction of larvae which
        #settle there.  note, this does mean there is some chance that
        #more will survivie then recruited their, but the average
        #should be right. Is this a problem?
        whichSpeciesSettle=searchsorted(thisCumsum,rng.random(Pmax))

        #breakpoint()
        
        #now put them back into P. Note there are two ways to do
        #this which should be benchmarked
        if False:
            P[nin,:]=bincount(whichSpeciesSettle,minlength=Nspecies)
        elif False:
            #P[nin,:]=0 #this is surprizingly slow, do elsewhere. 
            for n in whichSpeciesSettle:
                P[nin,n]+=1
        else:
            numbaFastCount(whichSpeciesSettle,P[nin,:])
                

        #breakpoint()
            
    return P

def runModelOnce(Pinit,R0,R1,Pmax,Tmax,nsp):
    '''

    Now the whole model is wrapped in this function. Pass in the
    intial size of the introduced species as Pinit, along with the
    larval production rate of the introduced (R0) and native(R1)
    species.

    Returns the number of generations run before the 0 population
    goes extinct (or we reach Tmax), and the final population distribution.

    nsp is the number of introduction we are running. It is passsed in only
    so we can print it out in diagnostics below.

    '''

    #record when starting for benchmarking
    ticBench=time.time()

    #how big is domain
    Ndomain=len(nxny2nlin)
    
    #==================================
    #now make initial P
    Nspecies=2
    P=zeros((Ndomain,Nspecies),dtype=int)
        
    #here we set the initial abundance of the native species as
    #Pmax-(abundance of introduced species) so that the initial
    #population density is at Pmax everywhere. 
    P[:,0]=Pinit
    P[:,1]=Pmax-P[:,0]
        
    #==================================
            
    #make lat and lon vectors for plotting
    lonVec=array([nlin2lonLat[n][0] for n in nlin2lonLat])
    latVec=array([nlin2lonLat[n][1] for n in nlin2lonLat])
        
    #now loop over time, and let species propogate
    #for now assume Nspecies=1
    #print('Starting time iterations with Ndomain size',Ndomain,'and Nspecies',Nspecies)
    for nt in range(Tmax):

        #plot before working, just to see initial condition
        #if true, plot
        if False:
            if remainder(nt,20)==0:
                clf()
                plot(lonVec,latVec,'k,')
                for nSp in range(Nspecies):
                    indx=P[:,nSp]>0
                    plot(lonVec[indx],latVec[indx],'.')
                #axis([-90,-50,24,60])
                title('Starting generation %d population %d'%(nt,sum(P.flatten())))
                draw()
                show(block=False)
                pause(0.1)
            
        tic=time.time()
        #now loop over each species, and figure out where its larvae will settle
        #serial solution
        whereSettle=findWhereSettle(P,0,Nspecies,R0,R1,Ndomain,Nspecies)

        #now there are two possibilities. Where the total number of
        #larvae reaching a location is less than or equal to Pmax,
        #they all survive. Where it is greater, choose Pmax survivors
        #randomly from the larvae. Should make sure number of recruits
        #in each species is not greater than number of larvae which
        #reach (really? or can the average number be right?). Do this
        #in a function so it is easy to profile. NOTE WELL, assume all
        #fitness differences are in the fecundity and dispersal, and
        #thus handled by findWhereSettle().
        #
        #calculate this with whichLarvaeSurvive(), which updates P
        tic2=time.time()
        Pold=P.copy()
        P=whichLarvaeSurvive(whereSettle,0,Ndomain,Pmax,Ndomain,Nspecies)
            

        now=time.time()
        if False:
            #print diagnostic every generation
            print('   generation %d took'%nt,now-tic,' and',
                  now-tic2,'; the population of the two species are',sum(P[:,0]),sum(P[:,1]),flush=True)

        #if the total population of 0 species (the introduced one) is 0, then quit
        if sum(P[:,0])==0:
            #print('species 0 whent extinct at nt=',nt)
            break

    #print what species and what process id -- only use when debugging to figure
    #out effciency issues. 
    #print('   done with region',nsp,'PID',os.getpid(),'in time',time.time()-ticBench,flush=True)
    
    #return how many generations the introduced species persisted, and what its final
    #distribution was
    return nt,P[:,0]
        
__spec__=None #ignore this line, it fixes a bug in multiprocessing on osx
if __name__=="__main__":


    #Because the results are stochastic, we need to run the model many
    #times.  The favored introduced species can only be introduced
    #into a single region at a time. So if Nregions is the number of
    #regions, and nRun is the number of times we want to introduce a
    #species into a single region in order to get good estimates of
    #the likelihood of persistence, than there will be Nregions*nRun
    #model runs.
    #
    #for speed, the Nregion model runs will be made in parallel. 
    nRun=100 ; print('FOR SPEED OF INITIAL USE, nRun IS SET TO 100. IT SHOULD BE LARGER IN MOST CASES')

    for nR in range(nRun):

        print('\nStarting run %d\n'%nR)
    
        #where find initial conditions, and fecundity of the two species
        initialConditionFile='initialConditions/'+ConnectivityModelName+'.zip'
        R0=16.0/2*1.1 #This is the number of larvae produced by the introduced species
        R1=16.0/2 #This is the number of larvae produced by the native species
        Pmax=1
        Tmax=601*4 

        #now, and important and subtle parameter. How many introductions
        #do we have? I.e. in each species location, how many adults do we
        #start with. Call this Nintro. If negative, fill entire species
        #range. Otherwise just release Nintro into the domain.
        Nintro=1

        #size of domain
        Ndomain=len(nxny2nlin)

        #now make the file name in which data will be saved. Then
        #check if the file exists, and if it does not exist, run the
        #model. This allows the code to interrupted and restarted.
        if Nintro<0:
            fileOut=('modelOutputRelativeFitness/twoSpecies_'+ConnectivityModelName
                     +'_Params_R0_%2.4f_R1_%2.4f_Tmax%3.3d_Pmax%2.2d_nRun%d.zip'%(R0,R1,Tmax,Pmax,nR))
        else:
            fileOut=('modelOutputRelativeFitness/twoSpecies_'+ConnectivityModelName
                     +'_Params_R0_%2.4f_R1_%2.4f_Tmax%3.3d_Pmax%2.2d_Nintro%2.2d_nRun%d.zip'%(R0,R1,Tmax,Pmax,Nintro,nR))

        #check if file exists, skips if it does not
        if os.path.isfile(fileOut):
            print('skipping because it exists:',fileOut)
        else:
            #run model

            if False:
                #simple debugging run
                def findClosest(lonVec,latVec,lonP,latP):
                    return argmin(((lonVec-lonP)/cos(deg2rad(latVec)))**2+(latVec-latP)**2)
                lonVec=array([nlin2lonLat[n][0] for n in nlin2lonLat])
                latVec=array([nlin2lonLat[n][1] for n in nlin2lonLat])
                n0=findClosest(lonVec,latVec,-74.42,39.21);
                Pinit=zeros((Ndomain,),dtype=int)
                Pinit[n0]=4

                Pfinal,nt=runModelOnce(Pinit,R0,R1,Pmax,Tmax)
            else:
                #load an initial condition file, and make a run for every discrete "species" in that
                #initial condition file.
                jnk=zarr.load(initialConditionFile)
                speciesList=jnk['nSpecies']
                Nregions=len(unique(speciesList))

                print('This run has',Nregions,'regions in which a novel species is introduced')
                print('note, progress bar only includes start of each species run...')

                ntVecFinal=zeros((Nregions,),dtype=int) #the number of generations each species lasts before going extinct
                finalPopVec=zeros((Ndomain,Nregions),dtype=int) #the distribution of each generation at the end of the run

                if False:
                    #serial run
                    for nsp in [0]: #range(Nregions): #assume species numbered in range(0,Nregions)
                        indx=speciesList==nsp
                        Pinit=zeros((Ndomain,),dtype=int)
                        Pinit[indx]=Pmax
                        assert False,'have not implimented Nintro yet'

                        nt,Pfinal=runModelOnce(Pinit,R0,R1,Pmax,Tmax,nsp)
                        ntVecFinal[nsp]=nt
                        finalPopVec[:,nsp]=Pfinal
                else:
                    #parallel run. Each introduction region will be run
                    #on a seperate core of the machine. 

                    #In general, in machines with
                    #hyperthreading, divide the number of CPU's
                    #reported by cpu_count() by two because that
                    #number includes virtual cores created by
                    #hyperthreading. On machines without
                    #hyperthreading, don't divide by two. But
                    #experiment...
                    nCPU=mp.cpu_count()#//2 

                    argVec=[]
                    for nsp in range(Nregions): #make a seperate run for each introduction region
                        indx=speciesList==nsp #find all points in species in initial range.

                        #Pinit is the number of introduced species at
                        #each location in the Ndomain points that make
                        #up the modeled habitat
                        Pinit=zeros((Ndomain,),dtype=int)

                        #now, figure out how many larvae to introduce
                        if Nintro<0:
                            #if Nintro<0, fill entire 
                            Pinit[indx]=Pmax
                        else:
                            #if Nintro>0, then fill the species range with the smaller
                            #of Nintro or Pmax*(number of locations in species initial range)
                            initialPoints=arange(len(speciesList))[indx] #indices of points in species initial range
                            #now make Pmax copies of these points; this is a list of all possible places an introduction
                            #could go
                            jnk=[]
                            for nP in range(Pmax):
                                jnk=jnk+list(initialPoints)
                            initialPoints=array(jnk)
                            shuffle(initialPoints) #do this to randomize where individuals are placed

                            #now randomly add points to the domain. We know we can't add too many
                            #because initialPoints's length is the maximum number of points an
                            #initial species range can hold.
                            numPoints=min(Nintro,len(initialPoints))
                            for thePoint in initialPoints[:numPoints]:
                                Pinit[thePoint]+=1

                            if Nintro>len(initialPoints):
                                print('For population',nsp,'Population initially saturated')

                        #argVec are the list of arguements to be fed into the multiprocessing routine
                        #which will run all the different introductions on different cores. 
                        argVec.append((Pinit,R0,R1,Pmax,Tmax,nsp))

                    tic=time.time()
                    with mp.Pool(nCPU) as pool:
                        #tqdm adds a progress bar (which is wonky, because it
                        #really tells you when function is called), but it is close. 
                        output=pool.starmap(runModelOnce,tqdm.tqdm(argVec,total=len(argVec))) #displays progress bar
                        #output=pool.starmap(runModelOnce,argVec) #run without progress bar
                    print('One run of introductions into all regions done in',time.time()-tic,'seconds',flush=True)

                    #save into ntVecFinal how long the introduced species persisted in region n
                    #save into finalPopVec what the final spatial distribution of the species introduced
                    #into region n was. 
                    for n in range(len(output)):
                        nt,Pfinal=output[n]
                        ntVecFinal[n]=nt
                        finalPopVec[:,n]=Pfinal


                #now save for each run what the results of introductions into each region was. 
                if True: #turn off for benchmarking
                    lonVec=array([nlin2lonLat[n][0] for n in nlin2lonLat])
                    latVec=array([nlin2lonLat[n][1] for n in nlin2lonLat])
                    zarr.save(fileOut,
                              lonVec=lonVec,latVec=latVec,finalPopVec=finalPopVec,ntVecFinal=ntVecFinal,speciesList=speciesList)


