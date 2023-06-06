import xarray as xr
import numpy as np
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f',dest='filename')
parser.add_argument('-s',dest='startYear',type=int)
parser.add_argument('-i',dest='initCond',default=None,help="file containing initial condition - will be added before the first time step of input file. If None, will duplicate first timestep of input file.")
parser.add_argument('-D',dest='outDir',default='members')
parser.add_argument('-a',dest='adjust',default=0,type=float,help="adjust time values by this number. Helpful e.g. for WACCM monthly means which are given at the end of the month.")
args = parser.parse_args()

filename = args.filename #'Dec5_F2000_HungaTungaEnsemble_monthly-0006.interp.nc'
startYear= args.startYear

ds = xr.open_dataset(filename,decode_times=False)
if 'member' in ds.dims:
    is_ens = True
    nmembs = len(ds.member)
else:
    is_ens = False
    nmembs = 1
    print('ASSUMING THERE IS NO MEMBER DIMENSION')

filename = filename.split('/')[-1]

if args.initCond is not None:
    initConds = xr.open_dataset(args.initCond,decode_times=False)


if not os.path.isdir(args.outDir):
    os.mkdir(args.outDir)

for y,yr in enumerate(range(startYear,startYear+nmembs)):
    if is_ens:
        memb = ds.isel(member=y)
    else:
        memb = ds
    units = memb.time.attrs['units']
    cal   = memb.time.attrs['calendar']
    # MiMA requires monthly timeseries, with 12 months/year, starting Jan
    #  thus, we add a full year to the start and the end
    if args.initCond is None:
        #time0 = memb.isel(time=0).assign_coords({'time':2*memb.time[0]-memb.time[1]})
        time0 = memb.isel(time=slice(0,12)).assign_coords({'time':memb.isel(time=slice(0,12)).time-360})
    else:
        time0 = xr.open_dataset(args.initCond)
    #timeE = memb.isel(time=-1).assign_coords({'time':2*memb.time[-1]-memb.time[-2]})
    timeE = memb.isel(time=slice(-12,None)).assign_coords({'time':memb.isel(time=slice(-12,None)).time+360})
    memb = xr.concat([time0,memb,timeE],'time')
    memb = memb.assign_coords({'time':memb.time+args.adjust})
    memb_year = units.split('since ')[1].split('-')[0]
    memb.time.attrs['units'] = units.replace(memb_year,'{:04d}'.format(yr))
    memb.time.attrs['calendar'] = '360_day'
    # now check that we don't have negative time values, as MiMA does not like that
    conv = xr.decode_cf(memb,use_cftime=True)
    #if conv.time.dt.year[0] < yr:
    #    nyr = yr-1
    #    memb = memb.assign_coords({'time':memb.time+360})
    #    memb.time.attrs['units'] = units.replace(memb_year,'{:04d}'.format(nyr))
    #memb.time.attrs['calendar'] = cal.replace('360_day','360')
    outFile = args.outDir+'/'+filename.replace('.nc','.{0:04d}.nc'.format(yr))
    if 'member' in conv:
        del conv['member']
    conv.to_netcdf(outFile,encoding={'time':{'units':'days since 0001-01-01 00:00:00','calendar':'360_day','dtype':'float'}},unlimited_dims=['time'])
    os.system("ncatted -a calendar,time,o,c,360 "+outFile)
    print(outFile)



