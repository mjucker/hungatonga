import xarray as xr
from aostools import climate as ac
from glob import glob
from dask.diagnostics import ProgressBar
import os


overwrite=False
parallel_pval=False

ctrldir = '/scratch/w40/mxj563/mima/aqua_sponge_5yr/'
pertdir = '/scratch/w40/mxj563/mima/aqua_sponge_pert_5yr/'


def AddMember(ds): 
        num = int(ds.encoding["source"].split('/')[-2]) 
        ds.coords['member'] = xr.DataArray([num],[('member',[num])],name='member') 
        ds = ds.assign_coords({'time':zero_time})
        return ds.resample(time='1M',label='left').mean()


def ReadMiMAEns(files): 
        ds = xr.open_mfdataset(files,decode_times=False,preprocess=AddMember,concat_dim='member')
        for var in ['time_bounds','average_T1','average_T2','average_DT']:
            if var in ds:
                del ds[var]
        if ds.time.attrs['calendar'] == '360': 
            ds.time.attrs['calendar'] = '360_day' 
        return xr.decode_cf(ds,use_cftime=True)

if overwrite:
        needfiles = {'ctrl':{
                'all' : True,
                'mean': True,
                'min' : True,
                'max' : True,
                'std' : True},
                     'pert': {
                'all' : True,
                'mean': True,
                'min' : True,
                'max' : True,
                'std' : True,
#                'pval': True,
                     }
                }
else:
        needfiles = {'ctrl':{'all':False},'pert':{'all':False}}
        basename = ctrldir.split('/')[-2]
        for meth in ['mean','max','min','std']:
                filename = basename+'_ens_{0}.nc'.format(meth)
                if os.path.isfile(filename):
                        needfiles['ctrl'][meth] = False
                else:
                        needfiles['ctrl'][meth] = True
                        needfiles['ctrl']['all']= True
        basename = pertdir.split('/')[-2]
        for meth in ['mean','max','min','std','pval']:
                filename = basename+'_ens_{0}.nc'.format(meth)
                if os.path.isfile(filename):
                        needfiles['pert'][meth] = False
                else:
                        needfiles['pert'][meth] = True
                        needfiles['pert']['all'] = True
                        if meth == 'pval':
                                needfiles['ctrl']['all'] = True


files_ctrl = glob(ctrldir+'00[3,4,5]?/*.atmos_davg.nc') 
files_pert = glob(pertdir+'00[3,4]?/*.atmos_davg.nc')

files_ctrl.sort()
files_pert.sort()

ens = {}

if needfiles['ctrl']['all']:
        zero_time = xr.open_dataset(ctrldir+'0030/0030.atmos_davg.nc',decode_times=False).time
        ens['ctrl'] = ReadMiMAEns(files_ctrl)
        ens['ctrl'].attrs['source'] = ctrldir.split('/')[-2]

if needfiles['pert']['all']:
        zero_time = xr.open_dataset(pertdir+'0030/0030.atmos_davg.nc',decode_times=False).time
        ens['pert'] = ReadMiMAEns(files_pert)
        ens['pert'].attrs['source'] = pertdir.split('/')[-2]

outFiles = {}
for key in ens.keys():
    outFiles[key] = {}
    basename = ens[key].attrs['source']
    for meth in needfiles[key].keys():
        if meth != 'all' and needfiles[key][meth]:
                outFiles[key][meth] = basename+'_ens_{0}.nc'.format(meth)
                if meth == 'mean':
                        delayed = ens[key].mean('member').to_netcdf(outFiles[key][meth],compute=False)
                elif meth == 'min':
                        delayed  = ens[key].min('member').to_netcdf(outFiles[key][meth],compute=False)
                elif meth == 'max':
                        delayed  = ens[key].max('member').to_netcdf(outFiles[key][meth],compute=False)
                elif meth == 'std':
                        delayed = ens[key].std('member').to_netcdf(outFiles[key][meth],compute=False)
                elif meth == 'pval':
                        pvals = []
                        for var in ens[key].data_vars:
                                ctrl_zm = ens['ctrl'][var].mean('lon').load()
                                pert_zm = ens['pert'][var].mean('lon').load()
                                pval = ac.StatTest(pert_zm,ctrl_zm,'KS',dim='member',parallel=parallel_pval)
                                pval.name = 'pval_'+var+'_zm'
                                pvals.append(pval)
                        delayed = xr.merge(pvals).to_netcdf(outFiles[key][meth],compute=False)
                with ProgressBar():
                        print(outFiles[key][meth])
                        delayed.compute()
    

