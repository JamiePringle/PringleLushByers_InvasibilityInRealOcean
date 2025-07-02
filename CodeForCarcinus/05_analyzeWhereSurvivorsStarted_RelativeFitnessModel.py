from pylab import *
from numpy import *
import zarr
import xarray as xr
import time
import multiprocessing as mp
import createLinearModel_module as cLM
import cartopy.crs as ccrs
import cartopy


#What model run are we using
ConnectivityModelName='E_CmaenasHab_depth1_minPLD40_maxPLD40_months5_to_6';
Pmax=1; Tmax=601*4; R0=8.8; R1=8.0; Nintro=1; nRun=1000

#load connectivity data
EfileName='transposes/twoSpecies_'+ConnectivityModelName+'.zarr'

#load initial conditions
initialConditionFile='initialConditions/'+ConnectivityModelName+'.zip'
jnk=zarr.load(initialConditionFile)
speciesList=jnk['nSpecies']

#load final answer
for nR in arange(nRun):
    if Nintro<=0:
        jnk=zarr.load('modelOutputRelativeFitness/twoSpecies_'+ConnectivityModelName
                      +'_Params_R0_%2.4f_R1_%2.4f_Tmax%3.3d_Pmax%2.2d_nRun%d.zip'
                      %(R0,R1,Tmax,Pmax,nR))
    else:
        jnk=zarr.load('modelOutputRelativeFitness/twoSpecies_'+ConnectivityModelName
                      +'_Params_R0_%2.4f_R1_%2.4f_Tmax%3.3d_Pmax%2.2d_Nintro%2.2d_nRun%d.zip'
                      %(R0,R1,Tmax,Pmax,Nintro,nR))

    print('reading run nR',nR,'number fixed',sum(jnk['ntVecFinal']==Tmax-1))
    if nR==0:
        #finalPopVec=jnk['finalPopVec']
        lonVec=jnk['lonVec']
        latVec=jnk['latVec']
        ntVecFinal=jnk['ntVecFinal']
        fixedVec=0*jnk['ntVecFinal']
        fixedVec[jnk['ntVecFinal']==Tmax-1]=1
    else:
        ntVecFinal=jnk['ntVecFinal']+ntVecFinal
        indx=jnk['ntVecFinal']==Tmax-1
        fixedVec[indx]=fixedVec[indx]+1

#normalize by nRun to get average ntVecFinal
ntVecFinal=ntVecFinal/nRun
fixedVec=fixedVec/nRun
        


#first plot where survivors started
Nspecies=2
Ndomain=len(lonVec)

figure(3,figsize=(7.98, 4.5)); clf(); style.use('ggplot')

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
ax1.add_feature(cartopy.feature.LAND,edgecolor='xkcd:light grey')
ax1.gridlines(draw_labels=True)

if True:
    #plot fraction fixed
    titleText='\nFraction persist for %d generations'%Tmax
    for n in range(len(unique(speciesList))):
        indx=speciesList==n
        if fixedVec[n]>0.5/100: #ignore less than some fraction
            scatter(lonVec[indx],latVec[indx],s=9,c=fixedVec[n]*100+0*lonVec[indx],vmin=0,
                    vmax=20.0+0*amax(fixedVec)*100,cmap='cool',zorder=round(fixedVec[n]*100+6),
                    transform=ccrs.PlateCarree())
        # else:
        #     if fixedVec[n]>0.0:
        #         #plot(lonVec[indx],latVec[indx],'k.',transform=ccrs.PlateCarree())
        #         scatter(lonVec[indx],latVec[indx],s=3,c='xkcd:grey',zorder=round(fixedVec[n]*100+6),
        #                 transform=ccrs.PlateCarree())
    

colorbar(shrink=0.7,location='left')
title('Nintro is %d, R0/R1=%4.2f'%(Nintro,R0/R1)+titleText,fontsize='small')

draw()
show(block=False)

