import xarray as xr
from dask.diagnostics import ProgressBar

file = 'waccm_pert_ens.nc'

ds = xr.open_dataset(file,chunks={})

da = ds.resample(time='QS-DEC').mean()

outFile = file.replace('waccm','waccm_season')
print(outFile)
delayed = da.to_netcdf(outFile,compute=False)
with ProgressBar():
    delayed.compute()

file = file.replace('pert','ctrl')
ds = xr.open_dataset(file,chunks={})

da0 = ds.resample(time='QS-DEC').mean()

outFile = file.replace('waccm','waccm_season')
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

