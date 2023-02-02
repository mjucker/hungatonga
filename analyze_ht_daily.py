import xarray as xr
from aostools import climate as ac
from glob import glob
from dask.diagnostics import ProgressBar
import os

nyears = 1
period = slice(None,360*nyears)

parallel_pval=True

ctrldir = '/g/data/w40/mxj563/work/HungaTonga/mima/aqua_sponge_10yr/'
pertdir = '/g/data/w40/mxj563/work/HungaTonga/mima/aqua_sponge_pert_5yr/'
#ctrldir = '/g/data/w40/mxj563/work/HungaTonga/mima/bench_SH/'
#pertdir = '/g/data/w40/mxj563/work/HungaTonga/mima/bench_SH_pert/'
pertdir = '/scratch/w40/mxj563/mima/aqua_sponge_pert_10yr/'


def AddMember(ds):
        if zero_time is not None:
                ds = ds.isel(time=period)
                num = int(ds.encoding["source"].split('/')[-2]) 
                ds.coords['member'] = xr.DataArray([num],[('member',[num])],name='member')
                ds = ds.assign_coords({'time':zero_time})
        return ds


def ReadMiMAEns(files): 
        ds = xr.open_mfdataset(files,decode_times=False,preprocess=AddMember,concat_dim='member')
        for var in ['time_bounds','average_T1','average_T2','average_DT']:
            if var in ds:
                del ds[var]
        if ds.time.attrs['calendar'] == '360': 
            ds.time.attrs['calendar'] = '360_day' 
        return xr.decode_cf(ds,use_cftime=True)

# daily output
files_ctrl = glob(ctrldir+'00[3,4,5]?/*.atmos_davg.nc')
files_ctrl.sort()
zero_time = None
ctrl_daily = ReadMiMAEns(files_ctrl)
ctrl_daily.attrs['source'] = ctrldir.split('/')[-2]

files_pert = glob(pertdir+'00[3,4,5]?/*.atmos_davg.nc')
files_pert.sort()
zero_time = xr.open_dataset(files_pert[0],decode_times=False).time
zero_time = zero_time.isel(time=period)
pert_daily = ReadMiMAEns(files_pert)
pert_daily.attrs['source'] = pertdir.split('/')[-2]


for ds in [ctrl_daily,pert_daily]:
        for var in ['zsurf','nv','latb','lonb']:
                del ds[var]

# convert ctrl into ensemble
ens = {}
ens['pert'] = pert_daily
ctrl_ens = []
for member in pert_daily.member:
        start = member.values+1
        stop  = member.values+nyears
        ctrl_tmp = ctrl_daily.sel(time=slice('{:04d}'.format(start),'{:04d}'.format(stop)))
        ctrl_tmp = ctrl_tmp.assign_coords({'time':pert_daily.sel(member=member).time})
        ctrl_tmp['member'] = member
        ctrl_ens.append(ctrl_tmp)
ens['ctrl'] = xr.concat(ctrl_ens,'member')
ens['ctrl'].attrs['source'] = ctrldir.split('/')[-2]


# write ensemble files

# full ensembles
for key in ens.keys(): 
         outFile = ens[key].attrs['source']+'_ens.nc'
         outFile = outFile.replace('pert','daily_pert')
         outFile = outFile.replace('ctrl','daily_ctrl')
         if not key in outFile:
                 outFile = outFile.replace('_ens.nc','daily_{0}_ens.nc'.format(key))
         delayed = ens[key].to_netcdf(outFile,compute=False) 
         print(outFile) 
         with ProgressBar(): 
             delayed.compute() 

