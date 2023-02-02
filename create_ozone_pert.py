import xarray as xr
import numpy as np

inFile = '/scratch/w40/mxj563/mima/bench_SH/INPUT/o3_monthly.nc'
o3_in = xr.open_dataset(inFile,decode_times=False)
nme = inFile.split('.nc')[-2].split('/')[-1]
latb = o3_in['latb']
lonb = o3_in['lonb']
phalf= o3_in['phalf']
o3_in.time.attrs['units'] = 'days since 0001-01-01'
o3_in.time.attrs['calendar'] = '360_day'
o3_in = xr.decode_cf(o3_in).o3_monthly

do3 = xr.open_dataarray('o3_diff_shifted.nc')
do3 = do3*48/29 
do3 = do3.rename(lev='pfull')
nyears = int(max(do3['time.year'].values))

do3int = do3.interp(pfull=o3_in.pfull,method='cubic',kwargs={'fill_value':0})
do3int = do3int.interp(lat=o3_in.lat,method='cubic',kwargs={'fill_value':'extrapolate'})

do3int = xr.concat([do3int.interp(time='0001-01-01',kwargs={'fill_value':0}),do3int],'time') 

if len(o3_in.lon) > 1:
    do3int = do3int.interp(lon=o3_in.lon,method='cubic',kwargs={'fill_value':'extrapolate'})
else:
    do3int = do3int.mean('lon').expand_dims({'lon':1})
do3int = do3int.transpose('time','pfull','lat','lon') 

o3_in = o3_in.assign_coords(time=do3int.time[:12])

def ShiftYears(ds,year_shift,daysperyear):
    from datetime import timedelta
    dt = ds.assign_coords(time=ds.time+timedelta(days=year_shift*daysperyear))
    dt.time.attrs = ds.time.attrs
    return dt

o3_all = [] 
for decade in [0,1,2]: 
    dt = ShiftYears(do3int,decade*nyears,365)
    o3_dec = []
    for year in np.unique(dt.time.dt.year): 
        dt_yr = dt.sel(time='{0:04d}'.format(year))
        o3_yr = dt_yr + ShiftYears(o3_in,int(year)-1,365)
        o3_dec.append(o3_yr)
    o3_all.append(xr.concat(o3_dec,'time'))
o3_comp = xr.concat(o3_all,'time')
enc = {'time': {'dtype': float,
                'units': 'days since 0001-01-01 00:00:00',
                'calendar': '360_day',
                }
       }
for coord in o3_comp.coords:
    o3_comp[coord].attrs = o3_in[coord].attrs
o3_comp.name = nme+'_delta'
o3_comp.attrs = o3_in.attrs
o3_comp = o3_comp.to_dataset()
o3_comp['lonb'] = lonb
o3_comp['latb'] = latb
o3_comp['phalf']= phalf
o3_comp.attrs = o3_in.attrs
o3_comp.attrs['method'] = 'ozone_pi plus waccm_F2000 hungatonga perturbation, repeated for each decade'
outFile = o3_comp.name+'.nc'
o3_comp.to_netcdf(outFile,encoding=enc,unlimited_dims=['time']) 
print(outFile)
print('YOU PROBABLY WANT TO CONVERT THE CALENDAR TO SOMETHING MIMA KNOWS:')
print("> ncatted -a calendar,time,o,c,'thirty_day_months' {}".format(outFile)) 
