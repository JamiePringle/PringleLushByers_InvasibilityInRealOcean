from pylab import *
from numpy import *
import zarr
import xarray as xr
import time
import os
import multiprocessing as mp
from multiprocessing import shared_memory
import createLinearModel_module as cLM

# This code models the stochastic population dynamics of many neutral
# species in a habitat whose connectivity was defined by the code in
# 00_makeConnectivityMatrices.py
#
# For speed, the conectivity matrices are converted from ones where
# location is defined by the (nx,ny) grid locations in the underlying
# circulation model to a linear domain defined by a single index. This
# is done for speed, and does not affect anything -- nothing in the
# underlying connectivity of the habitats is changed.

#Assumes the model output directory modelOutputNeutral exists
assert os.path.exists('modelOutputNeutral'), 'please create directory modelOutputNeutral to store model output'  

#Initialize random number generator
rng=np.random.default_rng()

#We load the large connectivity data in global scope, so we don't have
#to pass it into the parallel functions each time we call
#the. Sigh. The connectivity matrix specified by ConnectivityModelName
#was created by 00_makeConnectivityMatrices.py
print('Loading model',flush=True)
ConnectivityModelName='E_CmaenasHab_depth1_minPLD40_maxPLD40_months5_to_6'
EfileName='transposes/'+ConnectivityModelName+'.zarr'
EwhereTo,EnumTo,EfracReturn,Ecumsum,\
    nxny2nlin,nlin2nxny,nxny2lonLat,nlin2lonLat=cLM.makeLinearModel(EfileName)
print('       done loading model',flush=True)


#This function calculates how many larvae are produced by each species, and where
#they settle given the statistics in the connectivity matrices E.
def findWhereSettle(Pname,whereSettleName,lowerBound,upperBound,R,Ndomain,Nspecies):
    #see below for what these parameters are. Pname and whereSettleName are the
    #shared memory names for P and whereSettle

    tic=time.time()

    #The variables P, the population, and whereSettle, the number of
    #larvae which settle at each location for each species, are shared
    #memory arrays. Only processes those parts of the array that lie
    #between species lowerBound:upperBound,
    #e.g. P[:,lowerBound:upperBound] and only change
    #whereSettle[:,lowerBound:upperBound]. This is more complicated than
    #the usual mechanisms for sharing data while multiprocessing, but is
    #much, much faster in this case. 
    shm_P=shared_memory.SharedMemory(name=Pname)
    shm_whereSettle=shared_memory.SharedMemory(name=whereSettleName)
    P=ndarray((Ndomain,Nspecies),dtype=int,buffer=shm_P.buf)
    whereSettle=ndarray((Ndomain,Nspecies),dtype=int,buffer=shm_whereSettle.buf)

    #Only work on species that exist. To do this, find the total
    #population of each species by suming along the space axes. P has
    #dimensions [space,species]
    Psum=P.sum(axis=0)
    
    #now loop over species, and find where the larvae of each settle.
    for nsp in range(lowerBound,upperBound):
        if Psum[nsp]>0:
            #Pvec=P[:,nsp]

            #What is total number of larvae launched from each point which survive
            totalSurvive=R*P[:,nsp]*EfracReturn #note, can be much less than 1!

            #Now convert totalSurvive into an array of integers. The
            #likelyhood that it is rounded up is the fractional part,
            #the likelyhood that it is round down is the 1-fractional
            #part. I.e. 1.7 has a 30% chance of being rounded to 2,
            #and a 70% chance of being rounded to 1. We do this by
            #adding to totalSurvive a random uniform number between 0
            #and 1, and then rounding down.
            intSurvive=floor(rng.random(Ndomain)+totalSurvive).astype(int)

            #define thisWhereSettle, the number of larvae from the
            #species we are considering now in the loop which settle
            #at each location. So it is the length of the spatial
            #domain
            thisWhereSettle=zeros((Ndomain,),dtype=int)
            for n in arange(Ndomain)[P[:,nsp]>0]:
                #make array of integers which is list of which entries in EwhereTo of where the larvae are going
                whereGoList=searchsorted(Ecumsum[n],rng.random(intSurvive[n]))

                #now make array of what indices in the linear domain the larvae are going into
                whereSettleInDomain=EwhereTo[n][whereGoList]

                #whereSettle contains the number of larve from the species which reach
                #each point in the domain.
                if False:
                    #this is actually much slower, by a factor of two... But is this true
                    #for all model configurations? Perhaps experiment and see for your case
                    thisWhereSettle=thisWhereSettle+bincount(whereSettleInDomain,minlength=Ndomain)
                else:
                    for n in whereSettleInDomain:
                        thisWhereSettle[n]+=1


            #put answer into whereSettle
            whereSettle[:,nsp]=thisWhereSettle

    #print('      done with a processes in mp.Pool in',time.time()-tic,flush=True)
    return None #all communication done through shared memory arrays

#write a function that determines how many of the larvae in
#whereSettle survive to be adults in P.  Now there are two
#possibilities: Where the total number of larvae reaching a location
#is less than or equal to Pmax, they all survive. Where it is greater,
#choose Pmax survivors randomly from the larvae.  Do this in a
#function so it is easy to profile. NOTE WELL, assume all fitness
#differences are in the fecundity and dispersal, and thus handled by
#findWhereSettle() and not in whichLarvaeSurvive() 
def whichLarvaeSurvive(Pname,whereSettleName,lowerBound,upperBound,Pmax,Ndomain,Nspecies):
    #inputs:
    #  P: population matrix (Ndomain,Nspecies). Pass in name of shared memory in which it is stored
    #  whereSettle: (Ndomain,Nspecies) number of larvae of each species
    #     which reach each locations. Pass in name of shared memory in which it is stored
    #  Pmax, Ndomain, and Nspecies as defined below

    #NOTE, THIS CODE ASSUMES P HAS BEEN ZERO'D OUT BEFORE IT IS CALLED!!!!

    #the variables P, the population, and whereSettle, the number of
    #larvae which settle at each location, are shared memory
    #arrays. Only processes those parts of the array that lie between
    #species lowerBound:upperBound, e.g. P[:,lowerBound:upperBound] and only
    #change whereSettle[:,lowerBound:upperBound]
    shm_P=shared_memory.SharedMemory(name=Pname)
    shm_whereSettle=shared_memory.SharedMemory(name=whereSettleName)
    P=ndarray((Ndomain,Nspecies),dtype=int,buffer=shm_P.buf)
    whereSettle=ndarray((Ndomain,Nspecies),dtype=int,buffer=shm_whereSettle.buf)

    #find total number of recruits into each location, and put them in
    #population if the number is less or equal than Pmax. Not that all
    #arrays are of shape [spaceIndex,speciesIndex]
    totalRecruits=whereSettle.sum(axis=1)

    #find which locations this call to whichLarvaeSurvive() will deal with...
    #make an indx of all False's, then fill area we are interested in with True
    indxWorryAbout=zeros((Ndomain,),dtype=bool)
    indxWorryAbout[lowerBound:upperBound]=True
    
    #Fill all locations where the number of larvae attempting to settle in a region
    #is less than Pmax. This does NOT include case where no larvae get to region!
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

        #now put them back into P. Note there are two ways to do
        #this which should be benchmarked
        if False:
            P[nin,:]=bincount(whichSpeciesSettle,minlength=Nspecies)
        else:
            #P[nin,:]=0 #this is surprizingly slow, do elsewhere. 
            for n in whichSpeciesSettle:
                P[nin,n]+=1

    #return nothing, everything is communicated via shared memory
    return None

__spec__=None
if __name__=="__main__":

    #In order to better understand stochastic effects, lets run this
    #model nRun times, and save the results from each model run
    nRun=100; print('WARNING -- RUNNING ONLY 100 TIMES; FOR REAL SCIENCE RUN MORE TIMES')

    #first, loop over model runs -- because results are stochastic, you
    #want to run it multiple times. nRun controls how many times the model is run
    for nR in range(nRun):
        print('Starting run',nR)

        #Define initial population statistics. P array has size
        #(Ndomain,Nspecies), where Nspecies is the number of species
        #and Ndomain is the size of the domain. P is an integer array,
        #since it is the number of individuals in a habitat.
        #
        #The carrying capacity of the habitat is Pmax.
        #
        #Run model for Tmax generations.
        #
        #The number of larvae produced per adult which survives till
        #it could settle if in suitable habitat is R, which may be a
        #floating point number..
        Pmax=1
        Tmax=601
        R=16.0

        #define the intial condition file -- this is created by
        #01_makeInitalIntroductionRanges.py and delineates the regions
        #into which individual species are introduced.
        initialConditionFile='initialConditions/'+ConnectivityModelName+'.zip'

        #now, and important and subtle parameter. How many
        #introductions do we have? I.e. in each species location, how
        #many adults of a species we are following do we start
        #with. Call this Nintro. If negative, fill entire species
        #range defined by the initialConditionsFile. Otherwise just
        #release Nintro individuals into the domain, as long as each
        #region they are released into can hold more adults than
        #Nintro.  All other spaces in the domain are filled with
        #another species which we do not pay attention to. (It must be
        #filled with something, so that there is not an intial burst
        #of population growth to fill the empty habitat. This would
        #play with various statistics in funny ways.)
        Nintro=1
        
        #The arrays for the population P and where the larvae settle,
        #whereSettle, are created as shared memory arrays to reduce
        #overhead of communication when parallelization. Both must be
        #(Ndomain,Nspecies) in size, and thus cannot be made until we know
        #the number of species and the size of the domain.  so lets define
        #the size of the domain below (Ndomain), and then load the
        #initialization data and find the number of species. 

        #how big is domain
        Ndomain=len(nxny2nlin)

        #Now make initial P.  What is initial condition? The file
        #specified by initialConditionFile should have been created
        #with the same connectivity matrix as used in this model, and
        #is a vector of the same length as the linear domain, with a
        #series of integers from 0 to Nspecies-1 which indicate which
        #species start where. The initial condition is created by
        #01_makeInitialIntroductionRanges.py
        jnk=zarr.load(initialConditionFile)
        speciesList=jnk['nSpecies']

        #the number of species is the number of regions if Nintro<0,
        #and all regions start full. It is the number of regions+1 if
        #Nintro>0, because we need a species to fill otherwise empty
        #initial habitat. Otherwise the initial introduction number is
        #obscured by a burst of growth into otherwise empty habitats.
        #this filler species is species numbered
        #len(unique(speciesList))
        if Nintro<0:
            Nspecies=len(unique(speciesList))
        else:
            Nspecies=len(unique(speciesList))+1

        #how lets make the shared memory for the P and whereSettle arrays
        jnk=zeros((Ndomain,Nspecies),dtype=int) #this is only used to find size of shared memory to make

        shm_P=shared_memory.SharedMemory(create=True,size=jnk.nbytes)
        P=ndarray((Ndomain,Nspecies),dtype=int,buffer=shm_P.buf)
        P.fill(0)
        Pname=shm_P.name

        shm_whereSettle=shared_memory.SharedMemory(create=True,size=jnk.nbytes)
        whereSettle=ndarray((Ndomain,Nspecies),dtype=int,buffer=shm_whereSettle.buf)
        whereSettle.fill(0)
        whereSettleName=shm_whereSettle.name        

        #Lets make the initial distribution of species -- the initial
        #condition for this model run.  Loop over all the species, and
        #introduce them into the habitat at the appropriate location
        for nSp in unique(speciesList):
            indx=speciesList==nSp
            #P[indx,nSp]=Pmax #start habitat full

            #now, figure out how many larvae to introduce
            if Nintro<0:
                #if Nintro<0, fill entire 
                P[indx,NSp]=Pmax
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
                    P[thePoint,nSp]+=1

                #a quick sanity check
                if Nintro>len(initialPoints):
                    print('For population',nsp,'Population initially saturated')

                #now we don't want empty habitat in the domain. So any empty habitat will be occupied
                #by species Nspecies-1, which should be len(speciesList). So check all of the habitat patches
                #in this introduction region, and if any are not filled to Pmax, add that number of individuals
                #to nSpeciesOverflow
                nSpeciesOverflow=len(unique(speciesList))
                assert Nspecies==nSpeciesOverflow+1, 'oops, think about what silly thing you have done'
                initialPoints=arange(len(speciesList))[indx]
                for nHabPatch in initialPoints:
                    totalP=sum(P[nHabPatch,:-1])
                    if totalP<Pmax:
                        P[nHabPatch,nSpeciesOverflow]=Pmax-totalP
                    if totalP>Pmax:
                        assert False,'how did the total population at point exceed the carrying capacity?'

        #ok, now just double check that all habitat is filled to Pmax
        Ptotal=sum(P,axis=1)
        assert (Ptotal==Pmax).all(),'Not all of the initial habitat is at carrying capacity!?!?'

        #make lat and lon vectors for plotting
        lonVec=array([nlin2lonLat[n][0] for n in nlin2lonLat])
        latVec=array([nlin2lonLat[n][1] for n in nlin2lonLat])

        #Do we run in parallel?
        runParallel=True
        if runParallel:
            nCPU=mp.cpu_count()    #use for machines without hyperthreading (apple silicon)
            nCPU=mp.cpu_count()//2 #for machines with hyperthreading (most intel/amd machines)
            pool=mp.Pool(nCPU) #keep same pool for all runs, uncomment if using multiprocessing

            #do smaller chunks help parallelize more?
            #the parameter below is the number of jobs per
            #process in pool. 
            nChunksPerCPU=1 #no, increasing this does not help... 

        #now loop over time, and let species propogate
        #for now assume Nspecies=1
        print('Starting time iterations with Ndomain size',Ndomain,'and Nspecies',Nspecies)
        ticBench=time.time()
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
            #NOTE WELL. To make parallization more efficient, P and whereSettle are
            #shared with findWhereSettle() through shared memory arrays. So findWhereSettle()
            #has the silent side effect of changing whereSettle. It is important that findWhereSettle()
            #only changes the habitat it is responsible for. 
            whereSettle.fill(0)
            if not runParallel: #parallel or serial choice
                #serial solution
                findWhereSettle(Pname,whereSettleName,0,Nspecies,R,Ndomain,Nspecies)
            else:
                #parallel solution
                chunkVec=linspace(0,Nspecies,nCPU*nChunksPerCPU+1,dtype=int)
                argIn=[]
                for nsp in range(nCPU*nChunksPerCPU):
                    argIn.append((Pname,whereSettleName,chunkVec[nsp],chunkVec[nsp+1],R,Ndomain,Nspecies))

                #again, nothing in output, since findWhereSettle() changes whereSettle through
                #a shared memory array. 
                output=pool.starmap(findWhereSettle,argIn)

            if False:
                #stop after sometime and benchmark
                if nt==20:
                    print('after',nt,'iterations, total time is',time.time()-ticBench)
                    assert False,'stop now'

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
            P.fill(0) # zero out P, to avoid doing this slowly below
            if not runParallel:
                whichLarvaeSurvive(Pname,whereSettleName,0,Ndomain,Pmax,Ndomain,Nspecies)
            else:
                #parallel solution
                chunkVec=linspace(0,Ndomain,nCPU*nChunksPerCPU+1,dtype=int)
                argIn=[]
                for nsp in range(nCPU*nChunksPerCPU):
                    argIn.append((Pname,whereSettleName,chunkVec[nsp],chunkVec[nsp+1],Pmax,Ndomain,Nspecies))

                #again, nothing in output, since findWhereSettle() changes whereSettle through
                #a shared memory array. 
                output=pool.starmap(whichLarvaeSurvive,argIn)



            #if extinct, stop
            #if sum(whereSettle)==0:
            #    break

            #calculate how many species left
            numSpeciesLeft=sum(P.sum(axis=0)>0)

            #if the number of species is 1, check if only species left
            #is the filler species that is used when Nintro>0. If it
            #only the filler species is left, then bail from the run,
            #since it means that all of the introduced species have
            #gone extinct. Note that this means the software that
            #analyzes results must be able to cope with missing output
            #files.
            if Nintro>0:
                if numSpeciesLeft==1:
                    #only one species left. Check to see it is the filler species
                    numBySpecies=P.sum(axis=0)
                    if numBySpecies[-1]>0:
                        print('   Only filler species persists. No introduced species left. Bail from run')
                        break

            now=time.time()

            #now save model run every so often
            if remainder(nt,100)==0 and nt>0:
                #now make the file name in which data will be saved. 
                if Nintro<0:
                    fileOut=('modelOutputNeutral/manySpecies_'+ConnectivityModelName
                             +'_Params_R_%2d_Tmax%3.3d_Pmax%2.2d_nRun%d.zip'%(R,nt,Pmax,nR))
                else:
                    fileOut=('modelOutputNeutral/manySpecies_'+ConnectivityModelName
                             +'_Params_R_%2d_Tmax%3.3d_Pmax%2.2d_Nintro%2.2d_nRun%d.zip'%(R,nt,Pmax,Nintro,nR))
                print('   generation %d took'%nt,now-tic,' and',now-tic2,
                      'at species level, total species left',numSpeciesLeft)
                zarr.save(fileOut,
                          lonVec=lonVec,latVec=latVec,P=P)


        #close shared memory and multiprocessing pool
        if runParallel:
            pool.close()
            shm_P.close(); shm_P.unlink()
            shm_whereSettle.close(); shm_whereSettle.unlink()
