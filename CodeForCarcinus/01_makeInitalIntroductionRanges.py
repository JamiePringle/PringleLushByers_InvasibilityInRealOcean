from pylab import *
from numpy import *
import zarr
import time
import createLinearModel_module as cLM
from sklearn.neighbors import BallTree
from collections import Counter
import os

#the purpose of this code is to break the domain up into a larger
#number of small regions. The introductions of non-native species
#occurs in these regions. The reason we do this is to reduce the
#number of regions to examine -- it becomes computationally infeasable
#to study each habitat location seperately. However, if we make the
#regions too large, than it we will lose spatial resolution. In the
#code below, the size of the regions, in either kilometers or model
#habitat grid points, is defined. A certain amount of experimentation
#and oceanographic intuition will be helpful in choosing the optimal
#region size.

#if directory for the regions, "initialConditions/" does not exist, make it
if not os.path.exists('initialConditions'):
    os.mkdir('initialConditions')
    print('made directory "initialConditions" to store initial region locations')

#load the connectivity for the habitat you would like to process. 
ConnectivityModelName='E_CmaenasHab_depth1_minPLD40_maxPLD40_months5_to_6'
EfileName='transposes/'+ConnectivityModelName+'.zarr'
EwhereTo,EnumTo,EfracReturn,Ecumsum,\
    nxny2nlin,nlin2nxny,nxny2lonLat,nlin2lonLat=cLM.makeLinearModel(EfileName)
fileNameOut='initialConditions/'+ConnectivityModelName+'.zip'

#how big is domain
Ndomain=len(nxny2nlin)
    
#get domain locations in lat and lon
lonVec=array([nlin2lonLat[n][0] for n in nlin2lonLat])
latVec=array([nlin2lonLat[n][1] for n in nlin2lonLat])

#now make a plot that we will fill in
clf()
plot(lonVec,latVec,'k.')


#and here is where it gets arbitrary. For first try, pick a random
#point. find closest Ngroup points. Make one region. Iterate until
#no more points are not in a region.
if True:
    #this was what we did in points
    Nclose=5

#dictionary with key (nx,ny) to keep species number in. When done, len(nSpecies)==Ndomain
nSpecies={}


# ok, the algorythem is simple. heh.
#
# start with an initial habitat H which is a single point, and a
# dictionary of all points that are adjacent to it but not in any
# other species. Call this A
#
# START HERE: Add all points in A to to H, and make a newA with all points that
# are adjacent to points in A and not in H or any other species.
#
# make newA into A, and repeat, until size of H exceeds Nclose.
#
# then mark all the elements of H in nSpecies with a unique integer 
#
# now we need to pick a new intial point for H. First look in the last
# A, and see if there are any adjacent points not in any species. If
# yes, take the first one as the new point for a species area. If not,
# pick a random point not in any species.
#
# clear A and add the new point to a new initial H
#
# repeat until len(nSpecies)==Ndomain

#make a set of points in domain as (nx,ny)
habitatSet=set(nxny2nlin)

#make a set of all points that have not yet been assigned to a
#species. This starts as the entire habitat
notInSpecies=habitatSet.copy()

A=set()
a=notInSpecies.pop()
A.add(a) #start by adding something to A, the first point in the habitat.
notInSpecies.add(a) #it will be removed below, when the rest of A is added.
nsp=0 #what species are we working on. Will iterate upwards
while len(notInSpecies)>0:

    #print('starting species',nsp,'with',A)

    #this assumes that a an element of A
    lonLat=nxny2lonLat[a]
    
    H=set()
    while (len(H)<Nclose) and (len(A)>0):
        #Add A to H
        for a in A:
            H.add(a)

        #now find all neighbors of points in A and add to Anew if not already in H
        Anew=set()
        for a in A:
            #print('    working on',a)
            nx,ny=a
            for nnx in [-1,0,1]:
                for nny in [-1,0,1]:
                    #see if a point around a is not in H and is in habitat
                    if (not ((nx+nnx,ny+nny) in H)) and (((nx+nnx,ny+nny) in notInSpecies)) and ((nx+nnx,ny+nny) in habitatSet):
                        Anew.add((nx+nnx,ny+nny))
            #print('        done',flush=True)
        #print('   done with A',flush=True)
            
        #now move Anew to A, and repeat
        A=Anew

    #now add H to nSpecies, remove H from notInSpecies, and increase nsp by 1,
    #print('   starting to add H to species')
    for h in H:
        nSpecies[h]=nsp
        notInSpecies.remove(h)
    nsp+=1
    print('   added',len(H),'points to species',nsp,'with Nclose of',Nclose,flush=True)

    #ok, we need to pick a new point to start a new species. Lets pick one from Anew,
    #otherwise just get one from notInSpecies
    A=set()
    if len(Anew)>0:
        a=Anew.pop()
        A.add(a)        
        #print('Added new species from Anew',A,flush=True)
    else:
        if len(notInSpecies)>0:
            a=notInSpecies.pop()
            A.add(a) #start by adding something to A, the first point in the habitat.
            notInSpecies.add(a) #it will be removed below, when the rest of A is added.
            #print('Added new species from notInSpecies',A,flush=True)


#now make variable nSpeciesOut which is an array whose values give the
#species starting at each point in the linear habitat
nSpeciesOut=zeros((Ndomain,),dtype=int)
for k in nSpecies:
    nlin=nxny2nlin[k]
    nSpeciesOut[nlin]=nSpecies[k]

# OK, now for a reality check. There are a bunch of really small
# domains. Lets clump all less than some threshhold with adjacent
# species. When this is done, there will be no species with fewer than
# NcloseMin starting locations, if at all possible (which it might not
# be, if a species is on, e.g. a small island.
NcloseMin=2*Nclose//3
print('Now merging ranges which are too small. NcloseMin is',NcloseMin)
lenTooSmallRanges_old=-1 #see loop below, exit criterion
while True:
    #size of each species range in habitats
    speciesRangeSize=Counter(nSpeciesOut)
    tooSmallRanges=[k for k in speciesRangeSize if speciesRangeSize[k]<NcloseMin]
    tooSmallRangesSizes=[speciesRangeSize[k] for k in speciesRangeSize if speciesRangeSize[k]<NcloseMin]
    lenTooSmallRanges=len(tooSmallRanges)
    print('number of small ranges is now',lenTooSmallRanges,'number of ranges is',len(unique(nSpeciesOut)))

    if lenTooSmallRanges==0:
        print('stopped merging because we ran out of small species range sizes')
        break

    if lenTooSmallRanges==lenTooSmallRanges_old:
        print('stopped merging because we could not find any more species to merge')
        print(lenTooSmallRanges,'small ranges remain')
        break

    #just to make logic simpler, and because I don't care about speed
    #here, I will merge one species and then redo the loop. This
    #eliminates worrying about all sorts of edge cases.
    for speciesToMerge in tooSmallRanges:
        whereSpeciesToMerge=arange(Ndomain,dtype=int)[nSpeciesOut==speciesToMerge]

        #now go through all points adjacent to points in
        #whereSpeciesMerge and see if they are in a different
        #species. If so, merge and start again.
        pop2mergeTo=-1
        for aPoint in whereSpeciesToMerge:
            nx,ny=nlin2nxny[aPoint]
            for nnx in [-1,0,1]:
                for nny in [-1,0,1]:
                    if (nx+nnx,ny+nny) in nxny2nlin: #is in habitat
                        if nSpeciesOut[nxny2nlin[(nx+nnx,ny+nny)]]!=speciesToMerge:
                            pop2mergeTo=nSpeciesOut[nxny2nlin[(nx+nnx,ny+nny)]]
        if pop2mergeTo!=-1: #we found a different species!
            assert pop2mergeTo!=speciesToMerge,'oops, cant merge a species to itself!'
            nSpeciesOut[nSpeciesOut==speciesToMerge]=pop2mergeTo
            print('   merged',speciesToMerge,'with',pop2mergeTo,'producing range size of',sum(nSpeciesOut==pop2mergeTo))
            break #this break leaves "for speciesToMerge..." loop
        #else:
        #    print('   for',speciesToMerge,'found no neigbors with different species')
        #   

    #now remember what tooSmallRanges is
    lenTooSmallRanges_old=lenTooSmallRanges

#sigh, one last issue. the species must be number from 0..Nspecies-1, but after
#all the species merging, this is no longer true. So lets remap nSpeciesOut
#so that the species numbers go from 0 to Nspecies-1
speciesNumFuture=arange(len(unique(nSpeciesOut)),dtype=int)
nSpeciesOutFuture=0*nSpeciesOut
for nFuture,nNow in zip(speciesNumFuture,unique(nSpeciesOut)):
    indx=nSpeciesOut==nNow
    nSpeciesOutFuture[indx]=nFuture
assert len(unique(nSpeciesOutFuture))==len(unique(nSpeciesOut)),'number of species should be conserved'
assert amax(nSpeciesOutFuture)==len(unique(nSpeciesOutFuture))-1,'species should be numbered consequitively'
nSpeciesOut=nSpeciesOutFuture
    
#now plot
for n in unique(nSpeciesOut):
    indx=nSpeciesOut==n
    plot(lonVec[indx],latVec[indx],'.')
    draw()
    show(block=False)
    #pause(0.01)
    
#now save the species list
zarr.save(fileNameOut,nSpecies=nSpeciesOut)
