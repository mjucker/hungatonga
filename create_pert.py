import xarray as xr
import numpy as np

grid_file = 'mima/atmos_daily.nc'

restart_file = 'mima/spectral_dynamics.res.nc'


grid = xr.open_dataset(grid_file,decode_times=False)

# vertical profile from LMS
vertical_profile = [0,0.0899847,0.164989,0.211891,0.2258881,0.165305,0.0463565,0.00805761,0.0137417,0.0234326,0.0121541,0.00193966]
p1 = [100,82.5404,56.2341,38.3119,26.1016,17.7828,12.1153,8.25404,5.62341,3.83119,2.61016,1]
fprof = xr.DataArray(vertical_profile,coords=[('pfull',p1)],name='prof')
pprof = fprof.interp(pfull=grid.pfull).fillna(0)

restart = xr.open_dataset(restart_file)

## old version: 3D gaussian
#pres = grid.pfull
#z = -7.5*np.log(pres/1000)
##z = np.arange(len(pres))+1
lon = grid.lon
lat = grid.lat

phi0 = -20.5
lam0 = 184.6
#z0 = 26

dphi = 5
dlam = 10
#dz = 5.5
amp = 4.25e-4

#gauss = amp*np.exp(-(z-z0)**2/dz**2)*np.exp(-(lat-phi0)**2/dphi**2)*np.exp(-(lon-lam0)**2/dlam**2)

# new version: use MLS profile
gauss = amp*pprof*np.exp(-(lat-phi0)**2/dphi**2)*np.exp(-(lon-lam0)**2/dlam**2)

gauss.name = 'qa'
gaussp = gauss.to_dataset()
gaussp['lonb'] = grid['lonb']
gaussp['latb'] = grid['latb']
gaussp['phalf']= grid['phalf']
gaussp['lam0'] = float(lam0)
gaussp['phi0'] = float(phi0)
#gaussp['z0']   = float(z0)
gaussp['dphi'] = float(dphi)
gaussp['dlam'] = float(dlam)
#gaussp['dz']   = float(dz)
gaussp['amp']  = float(amp)
#gaussp.attrs['method'] = 'Gaussian specific humidity anomaly. Centers are lam0 [deg] in longitude, phi0 [deg] in latitude and z0 [km] in height. Widths are dlam in longitude, dlam in latitude and dz in height. Amplitude is amp [kg/kg].'
gaussp.attrs['method'] = 'Gaussian specific humidity anomaly. Centers are lam0 [deg] in longitude, phi0 [deg] in latitude. Widths are dlam in longitude, dlam in latitude. Amplitude is amp [kg/kg], and vertical profile derived from MLS.'
gaussp.to_netcdf('gauss.nc')

# compute total mass
from aostools import constants as ai
from aostools import climate as ac
mass = ac.GlobalMass(gauss)
#mass_y = ac.GlobalAvgXr(gauss,[-90,90])
#mass_x = np.deg2rad(mass_y.integrate('lon'))
#mass_p = mass_x.integrate('pfull')*100
#mass = ai.a0**2/ai.g*mass_p
print('Total added water mass = {0:f}Tg'.format(mass.values/1e9))

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

