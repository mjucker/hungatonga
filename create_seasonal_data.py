import xarray as xr
from dask.diagnostics import ProgressBar
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
args = parser.parse_args()
model = args.model

file = model+'_pert_ens.nc'

ds = xr.open_dataset(file,chunks={})

# label right to get DJF into the year of Jan-Feb, not Dec
da = ds.resample(time='QS-DEC',label='right').mean()

outFile = file.replace(model,model+'_season')
print(outFile)
delayed = da.to_netcdf(outFile,compute=False)
with ProgressBar():
    delayed.compute()

file = file.replace('pert','ctrl')
ds = xr.open_dataset(file,chunks={})

da0 = ds.resample(time='QS-DEC',label='right').mean()

outFile = file.replace(model,model+'_season')
print(outFile)
delayed = da0.to_netcdf(outFile,compute=False)
with ProgressBar():
    delayed.compute()

delta = da - da0

outFile = outFile.replace('ctrl','delta')
print(outFile)
delayed = delta.to_netcdf(outFile,compute=False)
with ProgressBar():
    delayed.compute()

