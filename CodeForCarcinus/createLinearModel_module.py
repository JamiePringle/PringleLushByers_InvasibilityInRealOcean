from pylab import *
from numpy import *
import zarr
import xarray as xr
from collections import defaultdict,Counter
import time
import getEZfateFromOSN


# This code takes a connectivity matrix and calculates the spread of a
# single species in a habitat where only a single individual can be in
# each patch.  R is with respect to the number of particles released,
# not the number that return (if I can).
#
# This is a testbed for fast model I am trying to develop, where the
# connectivity matrix is converted from one accessed by (nx,ny) model
# locations to habitat locations stores as in a linear array. This is
# done for speed. Details below.

#the spatial structure of the underlying hydrographic model is needed.
#the following downloads it if need be, and then sets maskFile to the
#file with this information. 
maskFile='EZfateData/EZfateFiles/ext-PSY4V3R1_mesh_zgr.nc'
maskFile=getEZfateFromOSN.getFileFromOSN(maskFile)

#This function returns the linear space module and other useful mappings
#when given the file name of a connectivity matrix as made by 00_makeTransposeMatrices.py
def makeLinearModel(EfileName):

    # what transpose file is going to give connectivity to the model?
    # It is called Enxny to make clear it is in nx,ny model grid space
    Enxny=zarr.open(EfileName,'r')

    #now lets make a map from (nx,ny) to a linear number of points 0 to
    #Npnts. The linear array will have one more point that that -- a point
    #that does not have any dispersal from it, so propagules that disperse
    #to it will effectively die. This will make fractional mortality
    #easier to implement later on. THE KEY ASSUMPTION IS that there are no
    #points in nxTo,nyTo that are not also in nxFrom,nyFrom.
    nxny2nlin={} #this maps from (nx,ny) model grid space to a linear space
    nxVec=Enxny.nxFrom[:] #load at once for speed
    nyVec=Enxny.nyFrom[:] #these can then be deleted if memory is an issue
    for n in range(len(Enxny.nxFrom)):
        nxny2nlin[nxVec[n],nyVec[n]]=n

    #now use this to find number of unique points, and make mapping from
    nlin2nxny={nxny2nlin[nkey]:nkey for nkey in nxny2nlin} #map from n to (nx,ny)

    #now a quick sanity check
    assert len(nxny2nlin)==len(nxVec),'there seem to be duplicate nxFrom,nyFrom points in E matrix?'

    #now find total size of domain, and make "kill bucket", a domain
    #location that does not connect to anywhere else
    Ndomain=len(nxVec)
    nlin2nxny[Ndomain]=(-1,-1) #an impossible location

    #now make a map from (nx,ny) and (nlin) to (lon,lat)
    gridData=xr.open_dataset(maskFile)
    nav_lon=gridData.nav_lon.data #again, load for speed
    nav_lat=gridData.nav_lat.data #could delete to get memory
    nxny2lonLat={}
    for nxny in nxny2nlin: #iterate over keys
        nxny2lonLat[nxny]=(nav_lon[nxny[1],nxny[0]],nav_lat[nxny[1],nxny[0]])
    nlin2lonLat={nxny2nlin[g]:nxny2lonLat[g] for g in nxny2nlin}

    #ok, now compile the Enxny matrix into the linear address space, in
    #the linear space. There will be no equivalent to (nxFrom,nyFrom),
    #since the nFrom will just be the index into the linear arrays. There
    #will be several arrays:
    #
    # EwhereTo[nFrom]=array of arrays of where to
    # EnumTo[nFrom]=array of arrays of number of propagules to
    # EfracReturn[nFrom]=an array of floats; the fraction of larvae launched that return
    # Ecumsum[nFrom]=array of the cumilative sum of numTo/sum(numTo), useful for picking a
    #                random point weighted by numTo
    #

    #lets calculate the fraction of releases that return to shore this
    #assumes that dictionaries were created from Enxny in row order
    EfracReturn=array([sum(n) for n in Enxny.numTo[:]])/Enxny.numLaunched[:]

    #there are sum zeros in numLaunched having to do with silliness in the Mercator grid.
    #so replace nan's with 0
    EfracReturn[isnan(EfracReturn)]=0.0

    #lets make EwhereTo. Load all data into memory for speed
    #if we end up short of memory, we can delete these
    nxTo=Enxny.nxTo[:]
    nyTo=Enxny.nyTo[:]
    numTo=Enxny.numTo[:]
    EwhereTo=empty(len(nxTo),dtype=object) 
    EnumTo=empty(len(nxTo),dtype=object)   
    Ecumsum=empty(len(nxTo),dtype=object)   
    for n in range(len(nxTo)):
        linTo=array([nxny2nlin[p] for p in zip(nxTo[n],nyTo[n])],dtype=int)
        EwhereTo[n]=linTo
        EnumTo[n]=numTo[n]
        Ecumsum[n]=cumsum(numTo[n])/sum(numTo[n])

    #now return what we have created
    return EwhereTo,EnumTo,EfracReturn,Ecumsum,nxny2nlin,nlin2nxny,nxny2lonLat,nlin2lonLat

if __name__=="__main__":

    EfileName='transposes/E_all_the_americas_depth20_minPLD22_maxPLD22_months4_to_6.zarr'
    EwhereTo,EnumTo,EfracReturn,Ecumsum, \
        nxny2nlin,nlin2nxny,nxny2lonLat,nlin2lonLat=makeLinearModel(EfileName)

    
