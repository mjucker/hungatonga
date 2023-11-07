#!env python
import xarray as xr
import numpy as np
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
from matplotlib import pyplot as plt
from matplotlib import ticker
sns.set_context('talk')
sns.set_style('whitegrid')
colrs = sns.color_palette()

# this has to be the same code as in create_pert.py!
grid_file = 'atmos_daily.nc'

restart_file = 'RESTART/spectral_dynamics.res.nc'


grid = xr.open_dataset(grid_file,decode_times=False)

# vertical profile from MLS
vertical_profile = [0,0.0899847,0.164989,0.211891,0.2258881,0.165305,0.0463565,0.00805761,0.0137417,0.0234326,0.0121541,0.00193966]
p1 = [100,82.5404,56.2341,38.3119,26.1016,17.7828,12.1153,8.25404,5.62341,3.83119,2.61016,1]
fprof = xr.DataArray(vertical_profile,coords=[('pfull',p1)],name='prof')
pprof = fprof.interp(pfull=grid.pfull).fillna(0)

fig,ax = plt.subplots()
pprof.sel(pfull=slice(100,None)).plot(y='pfull',ax=ax)

ax.set_ylim(1,100)
ac.LogPlot(ax)
ticks = ticker.FixedLocator(list(np.arange(1,10))+list(np.arange(10,100,10)))
ax.yaxis.set_minor_locator(ticks)
ax.yaxis.grid(True,which='minor')
sns.despine()
ax.set_xlim(0,0.25)

ax.set_title('vertical SWV anomaly profile')
ax.set_xlabel('fraction of total SWV anomaly []')

outFile = 'vertical_SWV_profile.pdf'
fig.savefig(outFile,bbox_inches='tight',transparent=True)
print(outFile)


