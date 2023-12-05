#import cf, cfplot as cfp
import numpy as np
import scipy as sp
import sys
import os
from netCDF4 import Dataset
import netCDF4 as nc
import matplotlib
import xarray as xr
import matplotlib.pyplot as plt
import warnings

import matplotlib.pyplot as plt
import matplotlib.path as mpath
# Quick plot to show the results
from cartopy import config
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import scipy.ndimage as ndimage
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import addcyclic
from cartopy.util import add_cyclic_point
import matplotlib.pylab as pl
from scipy import stats
import glob
from os.path import exists
import dask
dask.config.set(**{'array.slicing.split_large_chunks': True})

path='/gws/nopw/j04/realproj/users/adittus/'
var=sys.argv[1]
runid=sys.argv[2]
print(var)
#print(runid)

def prep_intermediate_files_realproj_4d(runid,var,plev,lev):
    files=glob.glob('/gws/nopw/j04/realproj/users/adittus/cdds/u-'+runid+'/cdds_data/CMIP6/CMIP/UKESM1-0-LL/piControl/r*i1p1f2/round-*/output/ap5/*/'+var+'/*nc')
    if (exists(path+'intermediate_files/stable/'+var+'/') == False):
        os.makedirs(path+'intermediate_files/stable/'+var+'/')
    if (exists(path+'intermediate_files/stable/'+var+'/'+var+'_'+plev+'_'+runid+'.nc') == False):
        for f in np.arange(0,len(files)):
            procvar = xr.open_dataset(files[f])[var][:,lev,:,:]
            procvar.to_netcdf(path+'intermediate_files/stable/'+var+'/'+var+'_'+plev+'_'+str(f)+'_'+runid+'.nc',engine='netcdf4')
        ua = xr.open_mfdataset(path+'intermediate_files/stable/'+var+'/'+var+'_'+plev+'_?_'+runid+'.*nc',combine='nested')[var]
        ua.to_netcdf(path+'intermediate_files/stable/'+var+'/'+var+'_'+plev+'_'+runid+'.nc',engine='netcdf4')
        for f in glob.glob(path+'intermediate_files/stable/'+var+'/'+var+'_'+plev+'_?_'+runid+'.nc'):
            os.remove(f)

def prep_intermediate_files_realproj_3d(runid,var):
    files=glob.glob('/gws/nopw/j04/realproj/users/adittus/cdds/u-'+runid+'/cdds_data/CMIP6/CMIP/UKESM1-0-LL/piControl/r*i1p1f2/round-*/output/ap5/*/'+var+'/*nc')
    print(files)
    if (exists(path+'intermediate_files/stable/'+var+'/') == False):
        os.makedirs(path+'intermediate_files/stable/'+var+'/')
    if (exists(path+'intermediate_files/stable/'+var+'/'+var+'_'+runid+'.nc') == False):
        for f in np.arange(0,len(files)):
            procvar = xr.open_dataset(files[f])[var][:,:,:]
            procvar.to_netcdf(path+'intermediate_files/stable/'+var+'/'+var+'_'+str(f)+'_'+runid+'.nc',engine='netcdf4')
        ua = xr.open_mfdataset(path+'intermediate_files/stable/'+var+'/'+var+'_?_'+runid+'.*nc',combine='nested')[var]
        ua.to_netcdf(path+'intermediate_files/stable/'+var+'/'+var+'_'+runid+'.nc',engine='netcdf4')
        for f in glob.glob(path+'intermediate_files/stable/'+var+'/'+var+'_?_'+runid+'.nc'):
            print('removing '+f)
            os.remove(f)
            
def prep_intermediate_gmst_realproj(runid):
    files=glob.glob('/gws/nopw/j04/realproj/users/adittus/cdds/u-'+runid+'/cdds_data/CMIP6/CMIP/UKESM1-0-LL/piControl/r*i1p1f2/round-2/output/ap5/Amon/tas/*nc')
    if (exists(path+'intermediate_files/stable/gmst/') == False):
        os.makedirs(path+'intermediate_files/stable/gmst/')
    if (exists(path+'intermediate_files/stable/gmst/tas_gmst_'+runid+'.nc') == False):
        for f in np.arange(0,len(files)):
            procvar = xr.open_dataset(files[f])['tas'][:,:,:].groupby('time.year').mean('time')
            procvar.to_netcdf(path+'intermediate_files/stable/gmst/tas_gmst_'+str(f)+'_'+runid+'.nc',engine='netcdf4')
        tas = xr.open_mfdataset(path+'intermediate_files/stable/gmst/tas_gmst_?_'+runid+'.*nc',combine='nested')['tas']
        weights = np.cos(np.deg2rad(tas.lat))
        tas_gmst=tas.weighted(weights).mean(("lon", "lat"))
        tas_gmst.to_netcdf(path+'intermediate_files/stable/gmst/tas_gmst_'+runid+'.nc',engine='netcdf4')
        for f in glob.glob(path+'intermediate_files/stable/gmst/tas_gmst_?_'+runid+'.*nc'):
            os.remove(f)

def prep_scenario_files_4d(runid,scenario,var,plev,lev,MIPtable):
    files=glob.glob('/badc/cmip6/data/CMIP6/*/MOHC/UKESM1-0-LL/'+scenario+'/'+runid+'/'+MIPtable+'/'+var+'/gn/latest/*nc')
    print(files)
    if files:
        if (exists(path+'intermediate_files/'+scenario+'/'+var+'/') == False):
            os.makedirs(path+'intermediate_files/'+scenario+'/'+var+'/')
        if (exists(path+'intermediate_files/'+scenario+'/'+var+'/'+var+'_'+plev+'_'+runid+'.nc') == False):
            for f in np.arange(0,len(files)):
                procvar = xr.open_dataset(files[f])[var][:,lev,:,:]
                print(procvar)
                procvar.to_netcdf(path+'intermediate_files/'+scenario+'/'+var+'/'+var+'_'+plev+'_'+str(f)+'_'+runid+'.nc',engine='netcdf4')
            ua = xr.open_mfdataset(path+'intermediate_files/'+scenario+'/'+var+'/'+var+'_'+plev+'_?_'+runid+'.*nc',combine='nested')[var]
            print(ua)
            ua.to_netcdf(path+'intermediate_files/'+scenario+'/'+var+'/'+var+'_'+plev+'_'+runid+'.nc',engine='netcdf4')
            for f in glob.glob(path+'intermediate_files/'+scenario+'/'+var+'/'+var+'_'+plev+'_?_'+runid+'.nc'):
                os.remove(f)

def prep_scenario_files_3d(runid,scenario,var):
    print(runid);print(scenario);print(var)
    files=glob.glob('/badc/cmip6/data/CMIP6/*/MOHC/UKESM1-0-LL/'+scenario+'/'+runid+'/*mon/'+var+'/g?/latest/*nc')
    print(files)
    if files:
        if (exists(path+'intermediate_files/'+scenario+'/'+var+'/') == False):
            os.makedirs(path+'intermediate_files/'+scenario+'/'+var+'/')
        if (exists(path+'intermediate_files/'+scenario+'/'+var+'/'+var+'_'+runid+'.nc') == False):
            for f in np.arange(0,len(files)):
                procvar = xr.open_dataset(files[f])[var][:,:,:]
                procvar.to_netcdf(path+'intermediate_files/'+scenario+'/'+var+'/'+var+'_'+str(f)+'_'+runid+'.nc',engine='netcdf4')
            if (scenario != 'piControl'):
                ua = xr.open_mfdataset(path+'intermediate_files/'+scenario+'/'+var+'/'+var+'_?_'+runid+'.*nc',combine='nested')[var]
                ua.to_netcdf(path+'intermediate_files/'+scenario+'/'+var+'/'+var+'_'+runid+'.nc',engine='netcdf4')
                for f in glob.glob(path+'intermediate_files/'+scenario+'/'+var+'/'+var+'_?_'+runid+'.nc'):
                    print('removing '+f)
                    os.remove(f)

def prep_scenario_gmst(runid,scenario):
    files=glob.glob('/badc/cmip6/data/CMIP6/*/MOHC/UKESM1-0-LL/'+scenario+'/'+runid+'/Amon/tas/gn/latest/*nc')
    if files:
        if (exists(path+'intermediate_files/'+scenario+'/gmst/') == False):
             os.makedirs(path+'intermediate_files/'+scenario+'/gmst/')
        if (exists(path+'intermediate_files/'+scenario+'/gmst/tas_gmst_'+runid+'.nc') == False):
            for f in np.arange(0,len(files)):
                var = xr.open_dataset(files[f])['tas'][:,:,:].groupby('time.year').mean('time')
                var.to_netcdf(path+'intermediate_files/'+scenario+'/gmst/tas_gmst_'+str(f)+'_'+runid+'.nc',engine='netcdf4')
            tas = xr.open_mfdataset(path+'intermediate_files/'+scenario+'/gmst/tas_gmst_?_'+runid+'.nc',combine='nested')['tas']
            weights = np.cos(np.deg2rad(tas.lat))
            tas_gmst=tas.weighted(weights).mean(("lon", "lat"))
            tas_gmst.to_netcdf(path+'intermediate_files/'+scenario+'/gmst/tas_gmst_'+runid+'.nc',engine='netcdf4')
            for f in glob.glob(path+'intermediate_files/'+scenario+'/gmst/tas_gmst_?_'+runid+'.nc'):
                os.remove(f)

# need to implement calculating "base" for all variables in function --> intermediate files - done
# need to implement seasonal analysis - done but check results are sensible
# need to implement capability for other variables - done
# consolidate functions for 3d and 4d vars
    
def get_seasindex_3(seas):
    season={'DJF':0,'JFM':1,'FMA':2,'MAM':3,'AMJ':4,'MJJ':5,'JJA':6,'JAS':7,'ASO':8,'SON':9,'OND':10,'NDJ':11,'ann':None}
    seasindex=season[seas]
    return(seasindex)

def calc_base(var,seas):
    seasindex=get_seasindex_3(seas)
    files=glob.glob(path+'intermediate_files/historical/'+var+'/'+var+'*_r*.nc')
    if files:
        if (seas == 'ann'):
            base=xr.open_mfdataset(files,combine='nested',concat_dim='cases')[var].groupby('time.year').mean('time')[:,0:50,:,].mean('cases').mean('year')
        else:
            rol=xr.open_mfdataset(files,combine='nested',concat_dim='cases')[var].rolling(min_periods=3, center=True, time=3).mean()[:,seasindex::12,:,:]
            base=rol[:,0:50,:,:].mean('cases').mean('time')
        return(base)

def find_indices_gwl_real(runid,gwl,var,seas,onesided_width=0.25):
    seasindex=get_seasindex_3(seas)
    files=glob.glob(path+'intermediate_files/stable/gmst/tas_gmst_'+runid+'.nc')
    if files:
        fileshist=glob.glob(path+'intermediate_files/historical/gmst/tas_gmst_r*.nc')
        filesreal=glob.glob(path+'intermediate_files/stable/'+var+'/'+var+'_*'+runid+'.nc')
#         print(files)
        gmst=xr.open_dataset(files[0])['tas'] # real gmst timeseries
        hist=xr.open_mfdataset(fileshist,combine='nested',concat_dim='cases')['tas'] # gmst timeseries for all hist runs --> baseperiod
        base=hist[:,0:50].mean('cases').mean('year')
        gmst=gmst-base
        mask=(gmst>=(gwl-onesided_width))&(gmst<=(gwl+onesided_width))
#         print(mask)
        if (seas == 'ann'):
            realdata=xr.open_dataset(filesreal[0])[var].groupby('time.year').mean('time') # read in realproj run to subset according to GWL
        else:
            rol=xr.open_dataset(filesreal[0])[var].rolling(min_periods=3, center=True, time=3).mean()
            realdata=rol[seasindex::12,:,:].groupby('time.year').mean('time')
        ua=realdata[mask,:,:]
        return(ua)

def find_indices_gwl_scenario(runid,gwl,var,seas,onesided_width=0.25):
    seasindex=get_seasindex_3(seas)
    files=glob.glob(path+'intermediate_files/ssp370/gmst/tas_gmst_'+runid+'.nc')
    if files:
        fileshist=glob.glob(path+'intermediate_files/historical/gmst/tas_gmst_r*.nc')
        filesreal=glob.glob(path+'intermediate_files/ssp370/'+var+'/'+var+'_*'+runid+'.nc')
        gmst=xr.open_dataset(files[0])['tas']
        hist=xr.open_mfdataset(fileshist,combine='nested',concat_dim='cases')['tas']
        base=hist[:,0:50].mean('cases').mean('year')
        gmst=gmst-base
        mask=(gmst>=(gwl-onesided_width))&(gmst<=(gwl+onesided_width))
        if (seas == 'ann'):
            realdata=xr.open_dataset(filesreal[0])[var].groupby('time.year').mean('time')
        else:
            rol=xr.open_dataset(filesreal[0])[var].rolling(min_periods=3, center=True, time=3).mean()
            realdata=rol[seasindex::12,:,:].groupby('time.year').mean('time')
#             realdata=realdata.rename({'time':'year'})
        ua=realdata[mask,:,:]
        return(ua)

def check_grid_for_var(var):
    filesreal=glob.glob(path+'intermediate_files/stable/'+var+'/'+var+'_*'+'bq777'+'.nc')
    lat=xr.open_dataset(filesreal[0])[var].lat
    lon=xr.open_dataset(filesreal[0])[var].lon
    return(len(lat),len(lon))
    
def extract_gwl_scenario(gwls,var,seas,return_means=True):
    latsize,lonsize=check_grid_for_var(var)
    print(latsize,lonsize)
    ua_scen=np.full((len(gwls),5000,latsize,lonsize),np.nan)
    ua_scen_means=np.full((len(gwls),latsize,lonsize),np.nan)
    for t in np.arange(0,len(gwls)):
        i=0
        for i in np.arange(1,21):
#             print('r'+str(i)+'i1p1f2')
            ua=find_indices_gwl_scenario('r'+str(i)+'i1p1f2',gwls[t],var,seas)
            if ua is not None:
                ua_scen[t,i:i+len(ua[:,0,0]),:,:]=ua
                i=i+len(ua[:,0,0])
        ua_scen_means[t,:,:] = np.nanmean(ua_scen[t,:,:,:],axis=0)
    if (return_means == True):
        return(ua_scen_means)
    else:
        return(ua_scen)

def extract_gwl_real(gwls,var,seas,return_means=True):
    latsize,lonsize=check_grid_for_var(var)
    runid = ('bq777','bw848','bz227','cd269','bw987','bu607')
    ua_real=np.full((len(gwls),5000,latsize,lonsize),np.nan)
    ua_real_means=np.full((len(gwls),latsize,lonsize),np.nan)
    for t in np.arange(0,len(gwls)):
        print('extracting GWL'+str(gwls[t]))
        i=0
        for run in runid: 
            #print(run)
            ua=find_indices_gwl_real(run,gwls[t],var,seas)
            ua_real[t,i:i+len(ua[:,0,0]),:,:]=ua
            i=i+len(ua[:,0,0])
        #del(i,ua)
        ua_real_means[t,:,:] = np.nanmean(ua_real[t,:,:,:],axis=0)
    if (return_means == True):
        return(ua_real_means,ua.lat,ua.lon)
    else:
        return(ua_real,ua.lat,ua.lon)

#runid = ('bq777','bw848','bz227','cd269','bw987','bu607')
for run in [runid]:
    print(run)
    #prep_intermediate_files_realproj_4d(run,var,plev='500',lev=5)
    #prep_intermediate_gmst_realproj(run)
    prep_intermediate_files_realproj_3d(run,var)
    
for i in [runid]: #np.arange(1,21):
    if ((var != 'ua') & (var != 'va') & (var != 'zg')):
        print('3d')
        #prep_scenario_files_3d('r'+str(i)+'i1p1f2','piControl',var)
        prep_scenario_files_3d('r'+str(i)+'i1p1f2','ssp245',var)
        # prep_scenario_files_3d('r'+str(i)+'i1p1f2','historical',var)
        # prep_scenario_files_3d('r'+str(i)+'i1p1f3','historical',var)
        # prep_scenario_gmst('r'+str(i)+'i1p1f2','ssp245')
        # prep_
#        scenario_gmst('r'+str(i)+'i1p1f3','ssp245')
    else:
        print('4d')
        prep_scenario_files_4d('r'+str(i)+'i1p1f2','ssp119',var,'850',5,'Amon')
        prep_scenario_files_4d('r'+str(i)+'i1p1f2','historical',var,'850',5,'Amon')
        prep_scenario_files_4d('r'+str(i)+'i1p1f3','historical',var,'850',5,'Amon')
