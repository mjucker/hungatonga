import xarray as xr
import numpy as np

grid_file = 'mima/atmos_daily.nc'

restart_file = 'mima/spectral_dynamics.res.nc'


grid = xr.open_dataset(grid_file,decode_times=False)


restart = xr.open_dataset(restart_file)

pres = grid.pfull
z = -7.5*np.log(pres/1000)
#z = np.arange(len(pres))+1
lon = grid.lon
lat = grid.lat

phi0 = -20.5
lam0 = 184.6
z0 = 28

dphi = 15
dlam = 15
dz = 3
amp = 2e-4

gauss = amp*np.exp(-(z-z0)**2/dz**2)*np.exp(-(lat-phi0)**2/dphi**2)*np.exp(-(lon-lam0)**2/dlam**2)

gauss.to_netcdf('gauss.nc')
# compute total mass
from aostools import constants as ai
from aostools import climate as ac
mass_y = ac.GlobalAvgXr(gauss,[-90,90])
mass_x = np.deg2rad(mass_y.integrate('lon'))
mass_p = mass_x.integrate('pfull')*100
mass = ai.a0**2/ai.g*mass_p
print('Total added water mass = {0:e}kg'.format(mass.values))

# put it into form for restart file
gauss = gauss.values
gauss = xr.DataArray(gauss,coords=[restart.zaxis_2,
                                          restart.yaxis_3,
                                          restart.xaxis_4],
                     name='sphum').expand_dims({'Time':2})

for attr,value in restart.sphum.attrs.items():
    gauss.attrs[attr] = value


restart_pert = restart.copy()

restart_pert['sphum'] = gauss
for attr,value in restart.sphum.attrs.items():
    restart_pert.sphum.attrs[attr] = value
for var in restart_pert.data_vars:
    if var != 'sphum':
        restart_pert[var].values = np.zeros_like(restart_pert[var])

restart_pert.to_netcdf('delta_init.nc')

