#!env python
import xarray as xr
import seaborn as sns
from aostools import climate as ac
from aostools import constants as at
from hungatonga import functions as fc
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from shapely import geometry
from cartopy import crs as ccrs
import numpy as np
import calendar
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
parser.add_argument('-v',dest='vars',default=None,nargs='+')
args = parser.parse_args()
model = args.model
sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()

if args.vars is None:
    fields = ['OLR','P','SLP','TS']
else:
    fields = args.vars

yrolls= [1,2,3,4,5]

var_fields = {'TS':'TREFHT'}

cmap_vars = {'TS':'t','P':'precip','OLR':'olr','SLP':'slp','default':'RdBu_r'}


scales = {'TS': 1, 'OLR': 1}

if 'waccm' in model:
    scales['P'] = 1e3*86400
    scales['SLP'] = 0.01
else:
    scales['P'] = 86400
    scales['SLP'] = 1

labls = {
    'TS': 'T2m [K]',
    'P' : 'Q [mm/day]',
    'OLR':'OLR [W/m2]',
    'SLP':'SLP [hPa]',
    }

# note: it is not practical to plot maps on the monthly rolling means
#        as this would create one figure per field per rolling mean with one panel per month
#        in the simulation, which is up to 120 months!

# create maps by season
sTS = xr.open_dataset('{0}_season_delta_ens.nc'.format(model))
sCT = xr.open_dataset('{0}_season_ctrl_ens.nc'.format(model))
#sTS = dTS.resample(time='QS-DEC').mean()
#sCT = CT.resample(time='QS-DEC').mean()

seasons = np.unique(sTS.time.dt.season)
# plot per variable per roll per season, 1 panel for each rolling year
for f,field in enumerate(fields):
    var = fc.variables[model][field]
    if var in var_fields.keys():
        var = var_fields[var]
    if field in cmap_vars.keys():
        cmap = at.cmaps[cmap_vars[field]]
    else:
        cmap = cmap_vars['default']
    da = sTS[var]*scales[field]
    dc = sCT[var]*scales[field]
    for a,roll in enumerate(yrolls):
        for s,seas in enumerate(seasons):
            filtr = sTS[var].time.dt.season == seas
            tmp = da.isel(time=filtr).rolling(time=roll).mean()
            ntimes = len(tmp.time.values)
            ncols = min(ntimes,3)
            nrows = (ntimes-1)//ncols+1
            fig,axs,transf = ac.Projection('PlateCarree',ncols=ncols,nrows=nrows,kw_args={'central_longitude':155})
            fig.set_figheight(nrows*3)
            fig.set_figwidth(ncols*4)
            std = dc.isel(time=filtr).rolling(time=roll).mean().std('member')
            tmp = tmp/std
            pval = ac.StatTest(tmp,0,'T','member',parallel=True)
            for t,time in enumerate(tmp.time):
                ax = axs.flatten()[t]
                tmp.mean('member').where(pval<0.1).isel(time=t).plot.contourf(levels=20,ax=ax,cmap=cmap,robust=True,**transf)
                ax.coastlines()
                ax.set_title('year {0}'.format(time.dt.year.values))
            fig.suptitle('{0}, {1}, {2}-year rolling mean'.format(labls[field],seas,roll))
            ac.AddPanelLabels(axs,'upper left',ypos=1.1) 
            outFile = 'figures/{0}_{1}_{2}_roll{3}_maps.pdf'.format(model,field,seas,roll)
            fig.savefig(outFile,bbox_inches='tight',transparent=True)
            print(outFile)
        plt.close('all')
