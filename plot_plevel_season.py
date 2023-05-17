import xarray as xr
from aostools import climate as ac
import seaborn as sns
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model')
parser.add_argument('-v',dest='var')
parser.add_argument('-l',dest='lev',type=int)
parser.add_argument('-L',dest='levels',default=None,type=float,nargs=3)
args = parser.parse_args()

seasons = ['DJF','JJA','MAM','SON']
olevel = args.lev
model = args.model
colrs = sns.color_palette()

var = args.var

#ctrls = xr.open_dataset(model+'_season_ctrl_ens.nc')[var].isel(time=slice(1,None))
diffs = xr.open_dataset(model+'_season_delta_ens.nc')[var].isel(time=slice(1,None)) #

dims = ac.FindCoordNames(diffs)
lev = dims['pres']
diffs = diffs.sel({lev:olevel})
#ctrls = ctrls.sel({lev:olevel})


mean = ['member','lat']

cwrap = 5
nrows = 10//5

if args.levels is None:
    levels = 20
else:
    levs = [*args.levels[:2],int(args.levels[2])]
    levs = [l for l in levs if l != 0 ]
    levels = np.linspace(*levs)

for s,season in enumerate(seasons):
    fig,axs,transf = ac.Projection('PlateCarree',ncols=cwrap,nrows=nrows,kw_args={'central_longitude':155})
    transf['add_colorbar'] = False
    fig.set_figheight(nrows*2*1.1)
    fig.set_figwidth(cwrap*4)
    filtr = diffs['time.season'] == season
    tmp = diffs.isel(time=filtr)
    pval = ac.StatTest(tmp,0,'T','member',parallel=True)
    for a,ax in enumerate(axs.flatten()):
        cf = tmp.where(pval<0.1).isel(time=a).mean('member').plot.contourf(ax=ax,x='lon',levels=levels,**transf)
        if a%cwrap > 0:
            ax.set_ylabel('')
        if a < cwrap:
            ax.set_xlabel('')
        ax.set_title('{3}{0}, {1}, year {2}'.format(olevel,season,a+1,var))
        ax.gridlines()
        ax.coastlines()
    ac.AddColorbar(fig,axs,cf)
    #fg.fig.suptitle(season)
    outFile = 'figures/{0}{1}_{2}.pdf'.format(var,olevel,season)
    fig.savefig(outFile,bbox_inches='tight')
    print(outFile)

    
