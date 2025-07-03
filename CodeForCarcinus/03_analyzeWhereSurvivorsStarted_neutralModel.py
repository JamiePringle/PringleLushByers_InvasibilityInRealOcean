from pylab import *
from numpy import *
import zarr
import xarray as xr
import time
import multiprocessing as mp
import createLinearModel_module as cLM
import cartopy.crs as ccrs
import cartopy
import os


#how many runs to include
nRun=100

#What model run are we using
ConnectivityModelName='E_CmaenasHab_depth1_minPLD40_maxPLD40_months5_to_6'; Pmax=1; Tmax=600; R=16.0; Nintro=1

#load connectivity data
EfileName='transposes/'+ConnectivityModelName+'.zarr'

#load initial conditions
initialConditionFile='initialConditions/'+ConnectivityModelName+'.zip'
jnk=zarr.load(initialConditionFile)
speciesList=jnk['nSpecies']


hasArunExisted=False
for nR in range(nRun):

    fileName=('modelOutputNeutral/manySpecies_'+ConnectivityModelName+
                  '_Params_R_%2.2d_Tmax%3.3d_Pmax%2.2d_Nintro%2.2d_nRun%d.zip'%(R,Tmax,Pmax,Nintro,nR))

    if not os.path.exists(fileName):
        print('does not exist:',fileName)
    else:
        #load final answer
        print('loading',nR)
        jnk=zarr.load(fileName)
        P=jnk['P']
        lonVec=jnk['lonVec']
        latVec=jnk['latVec']

        #what is P? It is the number of individuals of a given species at
        #a given location P[whichLocation,whichSpecies]. The total number
        #of individuals at each location, sum(P,axis=1) is no more than Pmax

        #calculate two things
        #
        #isPresent is the fraction of each species that exists anywhere and is Nspecies long
        #
        #Pave is [whichLocation,whichSpecies] long and gives the density
        #(population of species at point/Pmax) when averaged over all the
        #runs

        if not hasArunExisted:
            #for the first run we find, intialize variables
            hasArunExisted=True
            Pave=P.copy()
            totalPop=sum(P,axis=0)
            isPresent=zeros(shape(totalPop))
            indx=totalPop>0
            isPresent[indx]=isPresent[indx]+1
        else:
            Pave=Pave+P
            thisPop=sum(P,axis=0)
            totalPop=totalPop+thisPop
            indx=thisPop>0
            isPresent[indx]=isPresent[indx]+1

assert hasArunExisted,'Error, there was no ouput for any model run'
            
#normalize Pave and isPresent
Pave=Pave/nRun/Pmax
isPresent=isPresent
    
#====================================================================================
#first plot where survivors started
Nspecies=P.shape[1]
Ndomain=P.shape[0]

figure(1,figsize=(7.98, 7.5)); clf(); style.use('ggplot')

#define geographical range to display results over
lonMin=-80.0; lonMax=-48.0
latMin=35.0; latMax=55.0

#make map that we will use through whole run
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
    

#plot(lonVec,latVec,'r,',transform=ccrs.PlateCarree(),zorder=5)
for n in range(Nspecies):
    if totalPop[n]>0:
        indx=speciesList==n
        scatter(lonVec[indx],latVec[indx],s=15/2,c=isPresent[n]/nRun*100+0*latVec[indx],
                vmin=1.0,vmax=12.0,cmap='cool',transform=ccrs.PlateCarree(),zorder=round(isPresent[n]/nRun*100)+6)
colorbar(shrink=0.7,location='left')

title(EfileName+'\nstarting locations of survivors after %d generations'%Tmax,fontsize='small')
tight_layout()

draw()
show(block=False)

