import xarray as xr
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i',dest='inFile',help='Input file to shift time.')
parser.add_argument('-o',dest='outFile',help='Output file with shifted time.')
parser.add_argument('-s',dest='shift',type=int,help='Time shift in years.')
parser.add_argument('-d',dest='daysperyear',default=365,help='Number of days per year.')
args = parser.parse_args()

ds = xr.open_dataset(args.inFile,decode_times=False)
ds = ds.assign_coords(time=ds.time+args.shift*args.daysperyear)
ds.to_netcdf(outFile)


