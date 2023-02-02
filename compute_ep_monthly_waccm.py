import xarray as xr
from aostools import climate as ac
from aostools import constants as at
from glob import glob
import numpy as np
from dask.diagnostics import ProgressBar
from scipy.integrate import cumtrapz
import os

keep_vars = ['UV3d','VTH3d','UW3d','TH','U','T']

#ctrlfile = '/g/data/hh5/tmp/dd7103/Sept25-HungaTunga-PI/HungaTungaPI_monthly-06to36.nc'
#pertdir = '/g/data/hh5/tmp/dd7103/Hunga-Tunga-ensemble/merged-data/monthly/'
#member_length = 1

#ctrlfile = '/g/data/hh5/tmp/dd7103/F2000-HungaTungaEnsemble//control/F2000_HungaTungaControl_monthly-06to40.nc'
#pertdir = '/g/data/hh5/tmp/dd7103/F2000-HungaTungaEnsemble/'
#member_length = 5

ctrlfile = '/g/data/w40/mxj563/work/HungaTonga/waccm/F2000/F2000_HungaTunga_control_monthly-06to44yr.interp.nc'
pertdir  = '/g/data/w40/mxj563/work/HungaTonga/waccm/F2000/ensemble/'
#pertdir = '/scratch/w40/mxj563/hungatonga/tmp/'
member_length = 10
#ctrlfile = '/g/data/w40/mxj563/work/HungaTonga/waccm/PI/HungaTungaPI_monthly-06to36.interp.nc'
#pertdir  = '/g/data/w40/mxj563/work/HungaTonga/waccm/PI/ensemble/'
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

files = glob(pertdir+'*.nc')
files.sort()

tmp = xr.open_dataset(files[0],decode_times=False)
zero_time = DecodeDs(tmp).time

pert = xr.open_mfdataset(files,preprocess=AddMember,concat_dim='member')
pert = CleanVars(pert)
# restrict to member_length
pert = pert.isel(time=slice(None,12*member_length))

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


# compute EP fluxes
#  this code is borrowed from aostools
ep_comp  = {}
for key in ens.keys():
    if 'VTH3d' in ens[key].data_vars:
        vpThp = ens[key]['VTH3d'].mean('lon')
    else:
        factr = (1e3/ens[key].lev)**at.kappa
        vpThp = ens[key]['VT'].mean('lon')*factr
    filtr = vpThp<1e15
    vpThp = vpThp.where(filtr,0) #[m.K/s]
    if 'TH' in ens[key].data_vars:
        theta = ens[key]['TH'].mean(['member','lon'])
    else:
        theta = ens[key]['T'].mean(['member','lon'])*factr
    theta = theta.where(filtr,250) # [K]
    dTheta_dp = theta.differentiate('lev') # [K/hPa]
    #
    ubar = ens[key]['U'].mean('lon')
    fhat = np.deg2rad((ubar*at.coslat(ubar.lat)).differentiate('lat'))
    fhat = fhat/at.a0/at.coslat(ubar.lat)
    fhat = at.f(ubar.lat)-fhat # [1/s]
    #
    shear = ubar.differentiate('lev') # [m/s/hPa]
    #
    upvp = ens[key]['UV3d'].mean('lon')
    upvp = upvp.where(filtr,0) # [m2/s2]
    #
    ep1 = -upvp + shear*vpThp/dTheta_dp 
    ep1.name = 'ep1'
    ep1.attrs['units'] = 'm2/s2'
    #
    uw = ens[key]['UW3d'].mean('lon') # [m2/s2]
    uw = -uw.lev*at.g/at.Rd*uw/ens[key]['T'].mean('lon') # omega = - p.g/Rd*w/T
    uw = uw.where(filtr,0) # [m.hPa/s2]
    #
    ep2 = fhat*vpThp/dTheta_dp - uw
    ep2.name = 'ep2'
    ep2.attrs['units'] = 'hPa.m/s2'
    ep_comp[key] = xr.merge([ep1,ep2])
    #
    div1 = at.coslat(ep1.lat)*np.deg2rad(ep1.differentiate('lat')) - 2*at.sinlat(ep1.lat)*ep1
    div1 = div1/at.a0/at.coslat(ep1.lat)*86400
    div1.name = 'div1'
    div1.attrs['units'] = 'm/s/day'
    #
    div2 = ep2.differentiate('lev')*86400
    div2.name = 'div2'
    div2.attrs['units'] = 'm/s/day'
    #
    ep_comp[key] = xr.merge([ep1,ep2,div1,div2])
    

# write ensemble files
for key in ep_comp.keys():
    outFile = 'waccm_ep_{0}_ens.nc'.format(key)
    delayed = ep_comp[key].to_netcdf(outFile,compute=False)
    print(outFile)
    with ProgressBar():
        delayed.compute()

