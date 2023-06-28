import xarray as xr
from dask.diagnostics import ProgressBar
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
parser.add_argument('-v',dest='variables',default=None,nargs='+',help='Only work on this/these variable(s).')
args = parser.parse_args()
model = args.model

resArg = {'time':'QS-DEC','label':'left','loffset':'45D'}

def WorkFile(file):
    ds = xr.open_dataset(file,chunks={})
    if args.variables is not None:
        ds = xr.merge([ds[var] for var in args.variables])

    # label right to get DJF into the year of Jan-Feb, not Dec
    da = ds.resample(**resArg).mean()

    outFile = file.replace(model,model+'_season')
    if args.variables is not None:
        outStr = '_'.join(args.variables)
        outFile = outFile.replace('_season','_season_'+outStr)
    print(outFile)
    delayed = da.to_netcdf(outFile,compute=False)
    with ProgressBar():
        delayed.compute()
    return da,outFile

file = model+'_pert_ens.nc'
da,_ = WorkFile(file)

file = file.replace('pert','ctrl')
da0,outFile = WorkFile(file)

delta = da - da0

outFile = outFile.replace('ctrl','delta')
print(outFile)
delayed = delta.to_netcdf(outFile,compute=False)
with ProgressBar():
    delayed.compute()

