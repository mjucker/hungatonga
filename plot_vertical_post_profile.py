import xarray as xr
from hungatonga import functions as fc
from aostools import climate as ac
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib import ticker
import numpy as np

sns.set_style('whitegrid')
sns.set_context('talk')

lats = [35,50]
#lats = [-50,-35]
#years = [3,7]
years = '0001-10'
seasons = ['DJF','JJA']

ds = xr.open_dataset('waccm_season_delta_ens.nc',decode_times=False)
ds,_,_ = fc.CorrectTime(ds)

ds = ds.sel(lev=slice(None,500))
ds = ds['Q'].mean('lon')*1e6/0.622 #zonal mean in ppmv
ttle = 'zonal mean SWV anomalies,'
if isinstance(years,list):
    swv = ds.sel(time=slice('{0:04d}'.format(years[0]),'{0:04d}'.format(years[1]))).groupby('time.season').mean()
    swv = swv.sel(season=seasons)
    ssns = seasons
    ttle = ttle+'yrs {0}-{1}'.format(*years)
elif isinstance(years,str):
    swv = ds.sel(time=years).squeeze()
    ssns = None
    ttle = ttle+years

swv = ac.GlobalAvgXr(swv,lats)

pval = ac.StatTest(swv,0,'T','member')
swvm = swv.mean('member').where(pval<0.1)
std = swv.std('member').where(pval<0.1)

fig,ax = plt.subplots()
swvm.plot.line(y='lev',ax=ax)
colrs = sns.color_palette()
if ssns is not None:
    for s,seas in enumerate(swvm.season):
        ax.fill_betweenx(swvm.lev,(swvm-std).isel(season=s),(swvm+std).isel(season=s),color=colrs[s],alpha=0.3)
else:
    ax.fill_betweenx(swvm.lev,swvm-std,swvm+std,color=colrs[0],alpha=0.3)
#ax.set_xscale('log')
ax.set_ylim(ds.lev.values[0],1000)
ac.LogPlot(ax)
ticks = ticker.FixedLocator(list(np.arange(1,10))+list(np.arange(10,100,10)))
ax.yaxis.set_minor_locator(ticks)
ax.yaxis.grid(True,which='minor')
sns.despine()

if lats[0] < 0:
    lat0 = 'S'
else:
    lat0 = 'N'
if lats[1] < 0:
    lat1 = 'S'
else:
    lat1 = 'N'
ax.set_title(ttle+', {0}{2}-{1}{3}'.format(abs(lats[0]),abs(lats[1]),lat0,lat1))
ax.set_xlabel('water vapor anomaly [ppmv]')
