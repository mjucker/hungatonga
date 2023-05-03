import xarray as xr
from aostools import climate as ac
import seaborn as sns
import numpy as np

seasons = ['DJF','JJA','MAM','SON']
colrs = sns.color_palette()
cntrs = [l for l in np.linspace(-5,5,21) if l != 0]

ctrl = xr.open_dataset('waccm_ctrl_ens.nc')
pert = xr.open_dataset('waccm_pert_ens.nc')

ctrl = ctrl.FLNT
pert = pert.FLNT

diff = pert - ctrl

ctrls = ctrl.groupby('time.season').mean()
perts = pert.groupby('time.season').mean()

diffs = perts - ctrls

sel={'lat':slice(-20,20)}
mean = ['member']

fg = ctrls.sel(sel).mean(mean).plot(x='lon',col='season',col_wrap=2,cmap='RdBu')
fg.fig.set_figwidth(2*fg.fig.get_figwidth())

for a,ax in enumerate(fg.axs.flatten()):
    ttle = ax.get_title()
    season = ttle.split(' = ')[-1]
    tmp = diffs.sel(sel).sel(season=season)
    p = ac.StatTest(tmp,0,'T','member')
    (-tmp).mean(mean).where(p<0.1).plot.contour(levels=cntrs,x='lon',colors='black',ax=ax)
    ax.set_title(season)
    if a < 2:
        ax.set_xlabel('')
    if a%2 == 1:
        ax.set_ylabel('')
outFile = 'olr.pdf'
fg.fig.savefig(outFile,transparent=False,bbox_inches='tight')
print(outFile)


cwrap = 5

for s,season in enumerate(seasons):
    filtr = diff['time.season'] == season
    tmp = diff.isel(time=filtr).groupby('time.year').mean()
    tmp = tmp.sel(sel).mean('lat')
    fg=tmp.plot.line(x='lon',col='year',col_wrap=cwrap,color=colrs[s],alpha=0.2,add_legend=False)
    for a,ax in enumerate(fg.axs.flatten()):
        ax.axhline(0,color='k',lw=1,ls='--')
        tmp.isel(year=a).mean('member').plot(ax=ax,x='lon',add_legend=False,lw=2,color='k')
        if a%cwrap > 0:
            ax.set_ylabel('')
        if a < cwrap:
            ax.set_xlabel('')
        ax.set_title('{0}, year {1}'.format(season,a+1))
    fg.fig.suptitle(season)
    outFile = 'olr_{0}.pdf'.format(season)
    fg.fig.savefig(outFile,bbox_inches='tight')
    print(outFile)

    
