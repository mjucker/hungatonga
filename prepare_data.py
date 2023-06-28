###
# run this script to grab CAM4 data and bring it into form
###
import xarray as xr
import sys
import os
from glob import glob
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i',dest='inFile',help='The file containing data on hybrid grid.')
parser.add_argument('-o',dest='outFile',help='The file to write to.')
parser.add_argument('-O',dest='overwrite',action='store_true',help='Overwrite existing files.')
parser.add_argument('-t',dest='time_limit',default=None,type=int,nargs=2,help='Limit to given range of time steps. istart istop.')
args = parser.parse_args()

if not os.path.isfile('hybrid2pressure.template.ncl'):
    raise IOError('CANNOT PROCEED: NEED FILE hybrid2pressure.template.ncl IN WORKING DIRECTORY.')

if args.time_limit is not None:
    tsel = {'time':slice(args.time_limit[0],args.time_limit[1])}
else:
    tsel = None

overwrite = args.overwrite

dataFile = args.inFile
outFile3d  = args.outFile.replace('.nc','.3d.nc')
outFile2d  = args.outFile.replace('.nc','.2d.nc')


#vars_3d = ['O3','Z3','Q','U','V','T','RELHUM','UU','VV','VQ','VU','VT','VTH3d','QRS','QRL','OMEGA','CLOUD','CLDICE','CLDLIQ','TH']
vars_3d = ['O3','Z3','Q','U','V','T','TH','UV3d','VTH3d','UW3d','QRS','QRL','OMEGA','CLOUD','CLDICE','CLDLIQ']
#vars_3d = ['Q']
#vars_2d = ['CLDHGH','CLDLOW','CLDMED','CLDTOT','FLDS','FLDSC','FLNT','FLNTC','FSDS','FSDSC','FSNS','FSNSC','ICEFRAC','LANDFRAC','LHFLX','LWCF','OCNFRAC','PBLH','PHIS','PRECC','PRECL','PRECSC','PRECSL','PSL','SHFLX','SWCF','TMQ','TROPP_FD','TS','TSMN','TSMX','TREFHT',]
vars_2d = ['CLDTOT','CLDHGH','CLDLOW','CLDMED','FLNT','LANDFRAC','PHIS','PRECC','PRECL','PRECSC','PRECSL','PSL','TMQ','TROPP_FD','TS','TSMN','TSMX','TREFHT','U10','V10','FLDS','LWCF','SWCF','ICEFRAC']
#vars_2d = ['ICEFRAC']

## for MiMA input
#vars_3d = ['O3','Q','U','V','T','Z3']
#vars_2d = ['TS','PHIS']


def WriteFile(ds,name,enc=None):
    from dask.diagnostics import ProgressBar
    print('writing file',name)
    if 'time' in ds.coords and enc is None:
        enc={'time':{'dtype':float}}
    delayed = ds.to_netcdf(name,encoding=enc,compute=False)
    with ProgressBar():
        delayed.compute()

gridFile = 'level_grid.nc'

comp = None
atmos = None
if not os.path.isfile(outFile3d) or overwrite:
    comp = xr.open_dataset(dataFile,chunks={'time':1}).isel(tsel)
    dat3d = []
    for var in vars_3d:
        if var in comp and len(comp[var].dims) > 1:
            if 'lev' in comp[var].dims or 'ilev' in comp[var].dims:
                dat3d.append(comp[var])
    atmos = xr.merge(dat3d)
    WriteFile(atmos,outFile3d)

if not os.path.isfile(outFile2d) or overwrite:
    if comp is None:
        comp = xr.open_dataset(dataFile).isel(tsel)
    dat2d = []
    for var in vars_2d:
        if var in comp and len(comp[var].dims) > 1:
            dat2d.append(comp[var])
    SFC = xr.merge(dat2d)
    WriteFile(SFC,outFile2d)

if not os.path.isfile('data3.nc') or overwrite:
    if comp is None:
        comp = xr.open_dataset(dataFile).isel(tsel)
    out = xr.merge([comp.PS,comp.TS,comp.PHIS])
    WriteFile(out,'data3.nc')

if not os.path.isfile('level_grid.nc') or overwrite:
    if comp is None:
        comp = xr.open_dataset(dataFile).isel(tsel)
    grid = xr.merge([comp.hyam,comp.hybm,comp.P0,comp.hybi,comp.hyai])
    WriteFile(grid,'level_grid.nc')
print('DONE. YOU WILL NOW HAVE TO INTERPOLATE 3D DATA ONTO PRESSURE LEVELS.')
#######
# prepare bash script to do interpolation
##
# first, create ncl scripts
scripts = []
ncl_template_file = open('hybrid2pressure.template.ncl','r')
ncl_template = ncl_template_file.read()
if atmos is None:
    atmos = xr.open_dataset(outFile3d)
for var in atmos.data_vars:
            new_file = ncl_template.replace('INFILE',outFile3d)
            new_file = new_file.replace('OUTFILE','{0}.nc'.format(var))
            new_file = new_file.replace('VARNAME',var)
            if 'ilev' in atmos[var].coords:
                ya = 'hyai'
                yb = 'hybi'
            else:
                ya = 'hyam'
                yb = 'hybm'
            new_file = new_file.replace('HYA',ya).replace('HYB',yb)
            new_file = new_file.replace('GRIDFILE','level_grid.nc')
            new_file = new_file.replace('DATA3FILE','data3.nc')
            new_file = new_file.replace('LONGNAME',atmos[var].attrs['long_name'])
            new_file = new_file.replace('UNITS',atmos[var].attrs['units'])
            script_name = 'hybrid2pressure.{0}.ncl'.format(var)
            ncl_file = open(script_name,'w')
            ncl_file.write(new_file)
            ncl_file.close()
            scripts.append(script_name)
ncl_template_file.close()
# then, create bash script to call ncl scripts
bash_script = open('prepare_data.sh','w')
bash_script.write('#!/bin/bash\n')
bash_script.write('module load ncl nco\n')
bash_script.write('# interpolate onto pressure surfaces using NCL\n')
nscripts = len(scripts)
for s,script in enumerate(scripts):
    bash_script.write('echo "'+script+': {0}/{1}'.format(s+1,nscripts)+'"\n')
    bash_script.write('ncl -Q '+script+'\n')
vars =  list(atmos.data_vars)
bash_script.write('echo "COLLECT VARIABLES INTO ONE FILE"\n')
bash_script.write('mv {0}.nc {1}\n'.format(vars[0],args.outFile))
nvars = len(vars)-1
for f,field in enumerate(vars[1:]):
    bash_script.write('echo "{0}.nc: {1}/{2}"\n'.format(field,f+1,nvars))
    bash_script.write('ncks -A {0}.nc {1}\n'.format(field,args.outFile))
bash_script.write('echo "ADDING 2D FIELDS"\n')
bash_script.write('ncks -A {0} {1}\n'.format(outFile2d,args.outFile))
bash_script.write('echo "DONE. NOW YOU MIGHT WANT TO CLEAN UP TEMPORARY FILES WITH cleanup.sh"\n')
bash_script.close()
# finally, create a bash script to clean up intermediate files
clean_script = open('cleanup.sh','w')
clean_script.write('#!/bin/bash\n')
field_string = '{'+','.join(vars[1:])+'}'
clean_script.write('rm {0}.nc\n'.format(field_string))
clean_script.write('rm {0}\n'.format(outFile2d))
clean_script.write('rm {0}\n'.format(outFile3d))
clean_script.write('rm data3.nc\n'.format(outFile3d))
for script in scripts:
    clean_script.write('rm {0}\n'.format(script))
clean_script.write('rm level_grid.nc prepare_data.sh cleanup.sh\n'.format(outFile3d))
clean_script.write('echo "DONE."\n')
clean_script.close()
print('TO CONTINUE THE DATA PREPARATION, I HAVE NOW CREATED THE BASH SCRIPTS prepare_data.sh AND cleanup.sh')

