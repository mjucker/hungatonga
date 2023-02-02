import xarray as xr
from aostools import climate as ac
from glob import glob
from dask.diagnostics import ProgressBar
import os

#keep_vars = ['O3','U','V','T','TS','Q','PSL','TMQ','PRECC','PRECL','PRECSC','PRECSL','Z3','TCO','FLNT']
keep_vars = ['TS','TREFHT','Q','PSL','Z3','PRECC','PRECL','PRECSC','PRECSL','CLDICE']

sum_vars = {'PREC':['PRECC','PRECL','PRECSC','PRECSL']}


#ctrlfile = '/g/data/hh5/tmp/dd7103/Sept25-HungaTunga-PI/HungaTungaPI_monthly-06to36.nc'
#pertdir = '/g/data/hh5/tmp/dd7103/Hunga-Tunga-ensemble/merged-data/monthly/'
#member_length = 1

#ctrlfile = '/g/data/hh5/tmp/dd7103/F2000-HungaTungaEnsemble//control/F2000_HungaTungaControl_monthly-06to40.nc'
#pertdir = '/g/data/hh5/tmp/dd7103/F2000-HungaTungaEnsemble/'
#member_length = 5

ctrlfile = '/g/data/w40/mxj563/work/HungaTonga/waccm/F2000/F2000_HungaTunga_control_daily.interp.nc'
pertdir  = '/g/data/w40/mxj563/work/HungaTonga/waccm/F2000/daily/'
member_length = 30
#ctrlfile = '/g/data/w40/mxj563/work/HungaTonga/waccm/PI/HungaTungaPI_monthly-06to36.interp.nc'
#pertdir  = '/g/data/w40/mxj563/work/HungaTonga/waccm/PI/ensemble/'
#member_length = 2

def AddMember(ds):
    #num = int(ds.encoding["source"].split('.')[-2].split('-')[-1])
    num = int(ds.encoding["source"].split('.')[-3].split('-')[0])
    ds.coords['member'] = xr.DataArray([num],[('member',[num])],name='member')
    ds = ds.isel(time=slice(0,member_length))
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

# read daily output
ctrl = xr.open_dataset(ctrlfile)
ctrl = CleanVars(ctrl)
# add TCO
if 'O3' in keep_vars:
    ctrl['TCO'] = ac.TotalColumnOzone(ctrl['O3'],ctrl['T'])

files = glob(pertdir+'*.nc')
files.sort()

zero_time = xr.open_dataset(files[0]).isel(time=slice(0,member_length)).time
pert = xr.open_mfdataset(files,preprocess=AddMember,concat_dim='member')
pert = CleanVars(pert)
# restrict to member_length
pert = pert.isel(time=slice(None,12*member_length))
# add TCO
if 'O3' in keep_vars:
    pert['TCO'] = ac.TotalColumnOzone(pert['O3'],pert['T'])

if 'O3' in keep_vars:
    keep_vars.append('TCO')

ens = {'pert':pert}
# convert ctrl into ensemble
#  at the moment, we only have one member
#ctrl_ens = []
#for member in pert.member:
    #start = '{0:04d}-01-02'.format(member.values)
    #stop  = '{0:04d}-05-01'.format(member.values)
    #ctrl_tmp = ctrl.sel(time=slice(start,stop))
#    ctrl_tmp = ctrl.isel(time=slice(0,member_length))
#    if len(ctrl_tmp.time) == len(pert.sel(member=member).time):
#        ctrl_tmp = ctrl_tmp.assign_coords({'time':pert.sel(member=member).time})
#        ctrl_tmp['member'] = member
#        ctrl_ens.append(ctrl_tmp)
#    else:
#        print('CANNOT ADD MEMBER {0} AS TIME DIMENSION NOT SAME LENGTH: {1} ON CTRL, {2} ON PERT'.format(member.values,len(ctrl_tmp.time),len(pert.sel(member=member).time)))
#ens['ctrl'] = xr.concat(ctrl_ens,'member')
tmp = ctrl.isel(time=slice(0,member_length))
tmp = tmp.assign_coords({'time':pert.isel(member=0).time})
tmp['member'] = 0
ens['ctrl'] = tmp

# write ensemble files
for key in ens.keys():
    outFile = 'waccm_daily_{0}_ens.nc'.format(key)
    delayed = ens[key].to_netcdf(outFile,compute=False)
    print(outFile)
    with ProgressBar():
        delayed.compute()

