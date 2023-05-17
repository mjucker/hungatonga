import xarray as xr
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f',dest='filename')
parser.add_argument('-v',dest='var')
args = parser.parse_args()

filename = args.filename #'Dec5_F2000_HungaTungaEnsemble_monthly-0006.interp.nc'

var = args.var #= 'O3'

ds = xr.open_dataset(filename)[var]
ds.name = var.lower()
var = var.lower()
dsm= ds.convert_calendar('360_day',align_on='date')

# add edges
loni = list(dsm.lon.values)+[360]
lonbvals = np.array(loni)
lonvals = 0.5*(lonbvals[1:]+lonbvals[:-1])

pvals = dsm.lev.values
phvals = [0] + list(0.5*(pvals[1:]+pvals[:-1])) + [1000]

latvals = 0.5*(dsm.lat[1:].values+dsm.lat[:-1].values)
latb= dsm.rename({'lat':'latb'}).latb
latb.attrs = dsm.lat.attrs

lat = xr.DataArray(latvals,coords=[('lat',latvals)],name='lat')

lon = xr.DataArray(lonvals,coords=[('lon',lonvals)],name='lon')
lonb = xr.DataArray(lonbvals,coords=[('lonb',lonbvals)],name='lonb')
lonb.attrs = dsm.lon.attrs

dsm = dsm.rename({'lev':'pfull'})
dsm.pfull.attrs['units'] = 'hPa'
phalf = xr.DataArray(phvals,coords=[('phalf',phvals)],name='phalf')
phalf.attrs = dsm.pfull.attrs
dsm.pfull.attrs['edges'] = 'phalf'

ftype = 'f4'
ttype = ftype
enc = {'time':{'calendar':'360_day','dtype':ttype},
       'pfull':{'dtype':ftype},
       }

do = dsm.interp({'lon':lon,'lat':lat}).to_dataset()
do.lon.attrs = dsm.lon.attrs
do.lon.attrs['edges'] = 'lonb'
do.lat.attrs = dsm.lat.attrs
lat.attrs['edges'] = 'latb'

do['lonb'] = lonb
do['phalf']= phalf
do['latb'] = latb

do.encoding['unlimited_dims'] = ['time']

outFile = var+'.nc'
do.to_netcdf(outFile,encoding=enc)
print(outFile)
print('OBS: YOU WILL PROBABLY WANT TO CHANGE THE CALENDAR FROM 360_DAY TO 360 LIKE SO:')
print("ncatted -a calendar,{0},o,c,360 {1}".format(var,outFile)
