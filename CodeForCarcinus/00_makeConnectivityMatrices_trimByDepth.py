from pylab import *
from numpy import *
import zarr
import numcodecs
from collections import defaultdict,Counter
import makeConnectivityModule as mcm
import time
import popMatrixManipulate as pmm
import shutil
import os
import pickle
import bisect
import cartopy.crs as ccrs
import cartopy
import time
import multiprocessing as mp
import getEZfateFromOSN
import createLinearModel_module as cLM


#this function defines a habitat by water depth, as determined by
#the bathymetry on the Mercator GLORYS model grid.


#find where data is stored, in OSN data store
matDir='EZfateData/communityConnectivityMatrices/'

#make the function that will calculate the transpose, and take as
#argument what will be iterated over in the parallel runs. Look in the
#code to see what "inParam" does
def makeTranspose(inParam):

    #first, define the habitat that you will be examining. 

    #an easy way to define a polygon of lat lons is with
    #https://www.mapsdirections.info/en/draw-route-google-maps/
    #and save the KML file, and copy the lat/lon points from the KML file

    if True:
        #first, name your habitat -- this is your free choice
        habitatName='CmaenasHab'

        #what are the months in which dispersal happens?
        #for Carcinus maenas, in Maine, Most eggs are
        #released from May to June, Berrill 1982, and the
        #planktonic duration is about 40 days for mid-Maine
        #temperature -- Pringle et al. 2011
        inMonths=arange(5,6+1)
        inMonths=arange(5,6); print('#HACK FOR DEBUGING, ONLY LOAD ONE MONTH')

        #what depth are the larvae, in meters
        depth=1

        #what depths are in the habitat, in meters
        minHabDepth=0.0
        maxHabDepth=50.0

        #how long are they in the plankton; if inParam, the duration is
        #passed into the makeTranspose function
        minPLD=inParam; maxPLD=minPLD

        #what vertical behavior? "fixed" means fixed to a depth
        #"starts" allows the larvae to drift passively in all 3
        #dimensions
        vertBehavior='fixed'

        #make a polygon that encompasses your habitat. 
        regionPoly=[(-73.300781,34.537106),
                    (-79.013672,37.660994),
                    (-77.695313,41.32217),
                    (-73.740234,43.782001),
                    (-72.597656,46.990558),
                    (-67.412109,50.912558),
                    (-59.150391,52.438432),
                    (-51.328125,53.12947),
                    (-39.814453,50.690368),
                    (-44.296875,42.044194),
                    (-61.699219,35.257955)
                    ]
    else:
        assert False,'you need to choose one region'
        

    #define the name of the connectivity matrices, and where they will be stored, below
    #rootFileName defines the suffix of the filename. The connectivty forward in time
    #is in the file which starts with E, and the connectivity backwards in time is in
    #the matrices named Etranspose. 
    rootFileName='trimByDepth_%s_depth%d_minPLD%d_maxPLD%d_months%d_to_%d.zarr'%(habitatName,
                                                                      depth,minPLD,maxPLD,amin(inMonths),amax(inMonths))
    transposeFileName='transposes/Etranspose'+rootFileName
    connectFileName='transposes/E'+rootFileName

    #if the transpose file exists, do not re-create. If you do want to recreate it, delete the file
    if not os.path.exists(transposeFileName):
        #load a number of matrices to play with, and combine them to get a seasonal matrix
        store=zarr.MemoryStore()
        E=mcm.makeEmptyConnectivity(store)

        #loop over the regions that you wish to include -- there is no
        #need to download the data for a region if that region is not
        #in the polygon defined by regionPoly.
        for regionName in ['theAmericas']:#['AsiaPacific', 'EuropeAfricaMiddleEast', 'theAmericas']:
            for month in inMonths:
                print('   working on month',month,'in region',regionName,'for inParam',inParam,flush=True)
                matInFile=(matDir+
                           '%s/%dm/%s/'%(regionName,depth,vertBehavior)+
                           'climatology_month%2.2d_minPLD%2.2d_maxPLD%2.2d.zip'%(month,minPLD,maxPLD))

                #matInFile define locations of file that define the
                #model grid and Lagrangian connectivity from the
                #EZfate project, as described in
                #https://github.com/JamiePringle/EZfate the function
                #getEZfateFromOSN.getFileFromOSN() takes as an
                #argument a pathway to the EZfate data on the open
                #storage network S3 bucket, and returns the path to a
                #local file that contains the same data. Details as to
                #where and how this done, and where the data is
                #stored, can be found in the getEZfateFromOSN module.
                matInFile=getEZfateFromOSN.getFileFromOSN(matInFile)
                matIn=zarr.open(matInFile,'r')

                #if a habitat is defined by depth of habitat, then we
                #need a dictionary that gives the depth of the habitat at the
                #potential release depths
                gridDict=mcm.getDepthDict(matIn.nxFrom[:],matIn.nyFrom[:])
                
                #trim to all points within gridRadius of land
                #keep all true points, even if they fall outside of To, so we can
                #latter accurately count the number of points launched. 
                print('   trimming by distance from land')
                gridRadius=2.1*sqrt(2)
                matIn=mcm.trimConnectivity(matIn,gridDict,minHabDepth,maxHabDepth,keepAllTo=True)

                mcm.combineConnectivity(E,matIn)
        print('DONE reading in',inParam,'which has',E.nxFrom.shape[0],'points',flush=True)
        print(' ',flush=True)

        #now, add the numLaunched variable to E. This is the total number of To variables for each from variable
        numLaunched=[sum(p) for p in E.numTo]
        root=zarr.group(store=zarr.MemoryStore())
        numLaunchedZarrVar=root.empty(shape=(len(numLaunched),),name='numLaunched',dtype='i4')
        numLaunchedZarrVar[:]=array(numLaunched)
        zarr.convenience.copy(numLaunchedZarrVar,E,name='numLaunched')

        #breakpoint()
        
        #now trim to only include habitat inside of the regionPoly
        #polygon
        print('trimming by polygon')
        E=mcm.trimConnectivity_byPoly(E,regionPoly,keepInPoly=True)
        #breakpoint()

        #now trim connectivity again, but this time get rid of *To
        #points that land outside of habitat
        print('trimming To points')
        E=mcm.trimConnectivity(E,gridDict,minHabDepth,maxHabDepth,keepAllTo=False)
        
        #write out connectivity matrix
        #breakpoint()
        print('writing connectivity to disk',inParam,flush=True)
        if os.path.exists(connectFileName):
            shutil.rmtree(connectFileName)
        EoutStore=zarr.DirectoryStore(connectFileName)
        Eout=zarr.group(store=EoutStore)
        zarr.copy(E,Eout,name='/')
        print('   done writing',inParam,flush=True)

        #now invert matrix so we have the connectivity from where the
        #larvae end up to where they started (the inverse of what we
        #just calculated)
        EtransposeStore=zarr.MemoryStore() #in memory
        Etranspose=pmm.invertConMatrix(E,EtransposeStore)

        # normalizes the numTo array, which is the number of
        # lagrangian particles which go to a given point in E by the
        # sum of numTo, to get the likelyhood that the point in
        # (nxFrom,nyFrom) went to one the (nxTo,nyTo) points in the
        # matrix.  This normalized array is put in the variable
        # propTo.  Also calculates the cumilative sum of this matrix
        # and places it in cumSumPropTo

        print('normalizing inParam',inParam,flush=True)
        pmm.normalizeInvertConMatrix(Etranspose)
        print('   done normalizing',inParam,flush=True)

        #now write out to disk
        print('writing transpose to disk',inParam,flush=True)
        if os.path.exists(transposeFileName):
            shutil.rmtree(transposeFileName)
        EoutStore=zarr.DirectoryStore(transposeFileName)
        Eout=zarr.group(store=EoutStore)
        zarr.copy(Etranspose,Eout,name='/')
        print('   done writing',inParam,flush=True)
        print('\n\n')

    else:
        print('WARNING: skipping',transposeFileName,'because it exists already',flush=True)

    return connectFileName #return the fileName of the matrix we made, so we can plot it

    
#now run the function to make the habitat over a range of PLDs
for p in [40]:
    print('starting serial transpose for PLD of',p)
    connectFileName=makeTranspose(p)
    print('   done with',p)

    
#now, if the following is true, plot your habitat
if True:
    EwhereTo,EnumTo,EfracReturn,Ecumsum,\
        nxny2nlin,nlin2nxny,nxny2lonLat,nlin2lonLat=cLM.makeLinearModel(connectFileName)

    #how big is domain
    Ndomain=len(nxny2nlin)

    #get domain locations in lat and lon
    lonVec=array([nlin2lonLat[n][0] for n in nlin2lonLat])
    latVec=array([nlin2lonLat[n][1] for n in nlin2lonLat])

    #find the size of the domain, and then plot it
    lonMin=min(lonVec); lonMax=max(lonVec)
    latMin=min(latVec); latMax=max(latVec)

    #make map that we will use through whole run
    figure(1); clf()
    if True:
        #for limited area:
        central_lon=0.5*(lonMin+lonMax)
        central_lat=0.5*(latMin+latMax)
        proj=ccrs.Orthographic(central_lon, central_lat)

    ax1=subplot(1,1,1,projection=proj)
    ax1.set_extent((lonMin,lonMax,latMin,latMax))
    #ax1.coastlines()
    ax1.add_feature(cartopy.feature.LAND,edgecolor='Grey')
    ax1.gridlines(draw_labels=True)


    plot(lonVec,latVec,'r.',transform=ccrs.PlateCarree(),zorder=5,markersize=3)
    tight_layout()

    draw()
    show(block=False)
    
 
