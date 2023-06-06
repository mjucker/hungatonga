import xarray as xr
from aostools import climate as ac
from glob import glob
from dask.diagnostics import ProgressBar
import os

do_files = ['pert','ctrl']
#do_files = ['pert']

ctrldir = '/g/data/w40/mxj563/work/HungaTonga/mima/bench_SH_fix_ctrl/'
#pertdir = '/g/data/w40/mxj563/work/HungaTonga/mima/bench_SH_fix_125Tg/'
pertdir = '/g/data/w40/mxj563/work/HungaTonga/mima/bench_SH_fixq_o3init_125Tg/'
#pertdir = '/g/data/w40/mxj563/work/HungaTonga/mima/bench_SH_fixo3_qinit_125Tg/'

#ctrldir = '/g/data/w40/mxj563/work/HungaTonga/mima/bench_SH_ctrl/'
#pertdir = '/g/data/w40/mxj563/work/HungaTonga/mima/bench_SH_125Tg/'

#ctrldir = '/g/data/w40/mxj563/work/HungaTonga/mima/hthh_fix_ctrl/'
#pertdir = '/g/data/w40/mxj563/work/HungaTonga/mima/hthh_fix_125Tg/'

time_period = 10


def AddMember(ds): 
        if zero_time is not None:
                num = int(ds.encoding["source"].split('/')[-2]) 
                ds.coords['member'] = xr.DataArray([num],[('member',[num])],name='member')
                ds = ds.assign_coords({'time':zero_time})
        return ds


def ReadMiMAEns(files,addmem=True):
        if addmem:
                ds = xr.open_mfdataset(files,decode_times=False,preprocess=AddMember,combine='nested',concat_dim='member')
        else:
                ds = xr.open_mfdataset(files,decode_times=False) 
        for var in ['time_bounds','average_T1','average_T2','average_DT']:
            if var in ds:
                del ds[var]
        if ds.time.attrs['calendar'] == '360': 
            ds.time.attrs['calendar'] = '360_day' 
        return xr.decode_cf(ds,use_cftime=True)

# first the daily output
#files_ctrl = glob(ctrldir+'00[3,4,5,6]?/*.atmos_davg.nc')
#files_ctrl.sort()
#zero_time = None
#ctrl_daily = ReadMiMAEns(files_ctrl)

#files_pert = glob(pertdir+'00[3,4,5]?/*.atmos_davg.nc')
#files_pert.sort()
#zero_time = xr.open_dataset(files_pert[0],decode_times=False).time
#pert_daily = ReadMiMAEns(files_pert)

# then read monthly output
files_ctrl = glob(ctrldir+'00[3,4,5,6]?/*.atmos_monthly.nc')
files_ctrl.sort()
#zero_time = None
ctrl_monthly = ReadMiMAEns(files_ctrl,addmem=False)
#ctrl_monthly = xr.open_mfdataset(ctrldir+'00[3,4,5,6]?/*.atmos_monthly.nc',use_cftime=True)


files_pert = glob(pertdir+'00[3,4,5]?/*.atmos_monthly.nc')
files_pert.sort()
zero_time = xr.open_dataset(files_pert[0],decode_times=False).time
pert_monthly = ReadMiMAEns(files_pert)

## remove duplicate variables to make sure we can merge afterward
#bothvars = []
#for var in ctrl_monthly.data_vars:
#        if var in ctrl_daily.data_vars:
#                bothvars.append(var)
#for var in bothvars:
#        for ds in [ctrl_monthly,pert_monthly]:
#                del ds[var]
#for ds in [ctrl_monthly,ctrl_daily,pert_monthly,pert_daily]:
for ds in [ctrl_monthly,pert_monthly]:
        for var in ['zsurf','nv','latb','lonb']:
                del ds[var]


# now merge on monthly timescale
#ctrl_monthly = xr.merge([ctrl_monthly,ctrl_daily.resample(time='1M',loffset='16D',label='left').mean()])
#pert_monthly = xr.merge([pert_monthly,pert_daily.resample(time='1M',loffset='16D',label='left').mean()])
#pert_monthly = pert_monthly.transpose('member','time','pfull','lat','lon')
pert_monthly.attrs['source'] = pertdir.split('/')[-2]

# convert ctrl into ensemble
ens = {}
ens['pert'] = pert_monthly
ctrl_ens = []
for member in pert_monthly.member.values:
        start = member+1
        stop  = member+time_period
        ctrl_tmp = ctrl_monthly.sel(time=slice('{:04d}'.format(start),'{:04d}'.format(stop)))
        ctrl_tmp = ctrl_tmp.assign_coords({'time':pert_monthly.sel(member=member).time})
        ctrl_tmp['member'] = member
        ctrl_ens.append(ctrl_tmp)
ens['ctrl'] = xr.concat(ctrl_ens,'member')
ens['ctrl'].attrs['source'] = ctrldir.split('/')[-2]


# write ensemble files

# full ensembles
for key in ens.keys():
        if key not in do_files:
                continue
        outFile = ens[key].attrs['source']+'_ens.nc'
        if not key in outFile:
                outFile = outFile.replace('_ens.nc','_{0}_ens.nc'.format(key))
        delayed = ens[key].to_netcdf(outFile,compute=False) 
        print(outFile) 
        with ProgressBar(): 
                delayed.compute() 

