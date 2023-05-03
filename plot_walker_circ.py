import xarray as xr
from aostools import climate as ac
import seaborn as sns

seasons = ['DJF','JJA','MAM','SON']
olevel = 300
colrs = sns.color_palette()

ctrl = xr.open_dataset('waccm_ctrl_ens.nc')
pert = xr.open_dataset('waccm_pert_ens.nc')

ctrl = xr.merge([ctrl.U,ctrl.OMEGA*100])
pert = xr.merge([pert.U,pert.OMEGA*100])

ctrls = ctrl.groupby('time.season').mean()
perts = pert.groupby('time.season').mean()

diffs = perts - ctrls

sel={'lat':slice(-5,5),'lev':slice(100,None)}
mean = ['member','lat']

fg = ctrls.sel(sel).mean(mean).OMEGA.plot(x='lon',yincrease=False,col='season',col_wrap=2,cmap='RdBu')

for a,ax in enumerate(fg.axs.flatten()):
    ttle = ax.get_title()
    season = ttle.split(' = ')[-1]
    diffs.sel(sel).mean(mean).isel(lon=slice(None,None,5)).sel(season=season).plot.quiver(x='lon',y='lev',u='U',v='OMEGA',angles='xy',ax=ax,scale=5)
    ax.set_title(season)
    if a < 2:
        ax.set_xlabel('')
    if a%2 == 1:
        ax.set_ylabel('')
outFile = 'walker_circ.pdf'
fg.fig.savefig(outFile,transparent=False,bbox_inches='tight')
print(outFile)

do = pert.OMEGA - ctrl.OMEGA

cwrap = 5

for s,season in enumerate(seasons):
    filtr = do['time.season'] == season
    tmp = do.isel(time=filtr).groupby('time.year').mean()
    tmp = tmp.sel(lat=sel['lat'],lev=olevel).mean('lat')
    fg=tmp.plot.line(x='lon',col='year',col_wrap=cwrap,color=colrs[s],alpha=0.2,add_legend=False)
    for a,ax in enumerate(fg.axs.flatten()):
        ax.axhline(0,color='k',lw=1,ls='--')
        tmp.isel(year=a).mean('member').plot(ax=ax,x='lon',add_legend=False,lw=2,color='k')
        if a%cwrap > 0:
            ax.set_ylabel('')
        if a < cwrap:
            ax.set_xlabel('')
        ax.set_title('{0}hPa, {1}, year {2}'.format(olevel,season,a+1))
    #fg.fig.suptitle(season)
    outFile = 'omega{0}_{1}.pdf'.format(olevel,season)
    fg.fig.savefig(outFile,bbox_inches='tight')
    print(outFile)

    
