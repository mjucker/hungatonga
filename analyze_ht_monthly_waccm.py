import xarray as xr
from aostools import climate as ac
from glob import glob
from dask.diagnostics import ProgressBar
import os

keep_vars = ['O3','U','T','TS','Q','PSL','TMQ','PRECC','PRECL','Z3']

sum_vars = {'PREC':['PRECC','PRECL']}


ctrlfile = '/g/data/hh5/tmp/dd7103/Sept25-HungaTunga-PI/HungaTungaPI_monthly-06to36.nc'
pertdir = '/g/data/hh5/tmp/dd7103/Hunga-Tunga-ensemble/merged-data/monthly/'
member_length = 1

ctrlfile = '/g/data/hh5/tmp/dd7103/F2000-HungaTungaEnsemble//control/F2000_HungaTungaControl_monthly-06to40.nc'
pertdir = '/g/data/hh5/tmp/dd7103/F2000-HungaTungaEnsemble/'
member_length = 5

def AddMember(ds):
    num = int(ds.encoding["source"].split('.')[-2].split('-')[-1])
    ds.coords['member'] = xr.DataArray([num],[('member',[num])],name='member')
    ds = ds.isel(time=slice(None,-1))
    ds = ds.assign_coords({'time':zero_time})
    return ds

def CleanVars(ds):
    var_list = list(ds.data_vars)
    for var in var_list:
        if var in keep_vars:
            continue
        else:
            del ds[var]
    for var,sums in sum_vars.items():
        sumv = 0
        for sv in sums:
            sumv = sumv + ds[sv]
            del ds[sv]
        sumv.name = var
        ds = xr.merge([ds,sumv])
    return ds

# read monthly output
ctrl = xr.open_dataset(ctrlfile)
ctrl = CleanVars(ctrl)

files = glob(pertdir+'*.nc')
files.sort()

zero_time = xr.open_dataset(files[0]).isel(time=slice(None,-1)).time
pert = xr.open_mfdataset(files,preprocess=AddMember,concat_dim='member')
pert = CleanVars(pert)

ens = {'pert':pert}
# convert ctrl into ensemble
ctrl_ens = []
for member in pert.member:
    start = '{0:04d}-02'.format(member.values)
    stop  = '{0:04d}'.format(member.values+member_length-1)
    ctrl_tmp = ctrl.sel(time=slice(start,stop))
    if len(ctrl_tmp.time) == len(pert.sel(member=member).time):
        ctrl_tmp = ctrl_tmp.assign_coords({'time':pert.sel(member=member).time})
        ctrl_tmp['member'] = member
        ctrl_ens.append(ctrl_tmp)
    else:
        print('CANNOT ADD MEMBER {0} AS TIME DIMENSION NOT SAME LENGTH: {1} ON CTRL, {2} ON PERT'.format(member.values,len(ctrl_tmp.time),len(pert.sel(member=member).time)))
ens['ctrl'] = xr.concat(ctrl_ens,'member')

# write ensemble files
for key in ens.keys():
    outFile = 'waccm_{0}_ens.nc'.format(key)
    delayed = ens[key].to_netcdf(outFile,compute=False)
    print(outFile)
    with ProgressBar():
        delayed.compute()

