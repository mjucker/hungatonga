import xarray as xr
from aostools import climate as ac
import seaborn as sns

seasons = ['DJF','JJA','MAM','SON']
olevel = 500
colrs = sns.color_palette()

#ctrl = xr.open_dataset('waccm_ctrl_ens.nc')
#pert = xr.open_dataset('waccm_pert_ens.nc')

#ctrl = xr.merge([ctrl.U,ctrl.OMEGA*100])
#pert = xr.merge([pert.U,pert.OMEGA*100])

#ctrls = ctrl.groupby('time.season').mean()
#perts = pert.groupby('time.season').mean()

#diffs = perts - ctrls

ctrls = xr.open_dataset('waccm_season_ctrl_ens.nc').isel(time=slice(1,None))
diffs = xr.open_dataset('waccm_season_delta_ens.nc').isel(time=slice(1,None)) # remove first DJF
diffs = xr.merge([diffs.OMEGA*100,diffs.U])
ctrls = xr.merge([ctrls.OMEGA*100,ctrls.U])

sel={'lat':slice(-5,5),'lev':slice(100,None)}
mean = ['member','lat']

for season in seasons:
    seas = ctrls.time.dt.season == season
    fg = ctrls.isel(time=seas).sel(sel).mean(mean).OMEGA.plot(x='lon',yincrease=False,col='time',col_wrap=5,cmap='RdBu')

    for a,ax in enumerate(fg.axs.flatten()):
        #ttle = ax.get_title()
        #season = ttle.split(' = ')[-1]
        diffs.isel(time=seas).sel(sel).mean(mean).isel(time=a).OMEGA.plot.contour(x='lon',yincrease=False,colors='k',ax=ax)
        diffs.isel(time=seas).sel(sel).mean(mean).isel(lon=slice(None,None,5)).isel(time=a).plot.quiver(x='lon',y='lev',u='U',v='OMEGA',angles='xy',ax=ax,scale=5)
        ax.set_title('year {0}'.format(a+1))
        if a < 2:
            ax.set_xlabel('')
        if a%2 == 1:
            ax.set_ylabel('')
    outFile = 'figures/walker_circ_{0}.pdf'.format(season)
    fg.fig.savefig(outFile,transparent=False,bbox_inches='tight')
    print(outFile)

#do = pert.OMEGA - ctrl.OMEGA
do = diffs.OMEGA

cwrap = 5

for s,season in enumerate(seasons):
    filtr = do['time.season'] == season
    tmp = do.isel(time=filtr)#.groupby('time.year').mean()
    tmp = tmp.sel(lat=sel['lat'],lev=olevel).mean('lat')
    pval = ac.StatTest(tmp,0,'T','member',parallel=True)
    fg=tmp.plot.line(x='lon',col='time',col_wrap=cwrap,color=colrs[s],alpha=0.2,add_legend=False)
    for a,ax in enumerate(fg.axs.flatten()):
        ax.axhline(0,color='k',lw=1,ls='--')
        tmp.isel(time=a).mean('member').plot(ax=ax,x='lon',add_legend=False,lw=1,color='k',ls=':')
        tmp.where(pval<0.1).isel(time=a).mean('member').plot(ax=ax,x='lon',add_legend=False,lw=2,color='k')
        if a%cwrap > 0:
            ax.set_ylabel('')
        if a < cwrap:
            ax.set_xlabel('')
        ax.set_title('{0}hPa, {1}, year {2}'.format(olevel,season,a+1))
    #fg.fig.suptitle(season)
    outFile = 'figures/omega{0}_{1}.pdf'.format(olevel,season)
    fg.fig.savefig(outFile,bbox_inches='tight')
    print(outFile)

    
