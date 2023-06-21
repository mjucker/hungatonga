import xarray as xr
from aostools import climate as ac
from aostools import inout as ai
from glob import glob
from dask.diagnostics import ProgressBar
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-c',dest='case',default='100Tg')
args = parser.parse_args()

keep_vars = ['O3','U','V','T','TS','TREFHT','Q','PSL','PRECC','PRECL','PRECSC','PRECSL','Z3','TCO','OMEGA','FLNT','CLDICE','U10','V10','CLDTOT','FLDS','CLDHGH','CLDMED','CLDLOW','LWCF','SWCF']

sum_vars = {'PREC':['PRECC','PRECL','PRECSC','PRECSL']}

do_not_compress = ['Q','CLDICE']

case = args.case

sim_lengths = {'100Tg' : 10,
               '300Tg' : 5,
               '80Tg'  : 7
               }

#ctrlfile = '/g/data/hh5/tmp/dd7103/Sept25-HungaTunga-PI/HungaTungaPI_monthly-06to36.nc'
#pertdir = '/g/data/hh5/tmp/dd7103/Hunga-Tunga-ensemble/merged-data/monthly/'
#member_length = 1

#ctrlfile = '/g/data/hh5/tmp/dd7103/F2000-HungaTungaEnsemble//control/F2000_HungaTungaControl_monthly-06to40.nc'
#pertdir = '/g/data/hh5/tmp/dd7103/F2000-HungaTungaEnsemble/'
#member_length = 5

ctrlfile = '/g/data/w40/mxj563/work/HungaTonga/waccm/F2000/F2000_HungaTunga_control_monthly-06to44yr.interp.nc'
#ctrlfile = '/scratch/w40/mxj563/sh_met_book/F2000_HungaTonga_control_monthly-06to44yr.interp.nc'
pertdir  = '/g/data/w40/mxj563/work/HungaTonga/waccm/F2000/{0}/ensemble/'.format(case)
#pertdir  = '/g/data/w40/mxj563/work/HungaTonga/waccm/F2000/ensemble/high_top/'
#pertdir = '/scratch/w40/mxj563/hungatonga/tmp/'

member_length = sim_lengths[case]
#ctrlfile = '/g/data/w40/mxj563/work/HungaTonga/waccm/PI/HungaTungaPI_monthly-06to36.interp.nc'
#pertdir q = '/g/data/w40/mxj563/work/HungaTonga/waccm/PI/ensemble/'
#member_length = 2

def AddMember(ds):
    #num = int(ds.encoding["source"].split('.')[-2].split('-')[-1])
    num = int(ds.encoding["source"].split('.')[-3].split('-')[-1])
    ds.coords['member'] = xr.DataArray([num],[('member',[num])],name='member')
    #ds = ds.isel(time=slice(None,-1))
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

def DecodeDs(ds):
    attrs = ds.time.attrs
    ds = ds.assign_coords({'time':ds.time-15})
    ds.time.attrs = attrs
    return xr.decode_cf(ds)

# read monthly output
ctrl = xr.open_dataset(ctrlfile,decode_times=False)
ctrl = CleanVars(ctrl)
ctrl = DecodeDs(ctrl)

# add TCO
ctrl['TCO'] = ac.TotalColumnOzone(ctrl['O3'],ctrl['T'])

files = glob(pertdir+'*.nc')
files.sort()

tmp = xr.open_dataset(files[0],decode_times=False)
zero_time = DecodeDs(tmp).time

pert = xr.open_mfdataset(files,preprocess=AddMember,combine='nested',concat_dim='member')
pert = CleanVars(pert)
# restrict to member_length
pert = pert.isel(time=slice(None,12*member_length))
# add TCO
pert['TCO'] = ac.TotalColumnOzone(pert['O3'],pert['T'])

if 'TCO' not in keep_vars:
    keep_vars.append('TCO')

ens = {'pert':pert}
# convert ctrl into ensemble
ctrl_ens = []
for member in pert.member:
    start = '{0:04d}'.format(member.values)
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
    cvars = [var for var in ens[key].data_vars if var not in do_not_compress]
    enc = ai.DefCompress(ens[key],cvars)
    delayed = ens[key].to_netcdf(outFile,encoding=enc,compute=False)
    print(outFile)
    with ProgressBar():
        delayed.compute()

