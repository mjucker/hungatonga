import xarray as xr
from aostools import climate as ac
from hungatonga import functions as fc
import seaborn as sns
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model')
parser.add_argument('-v',dest='var')
parser.add_argument('-l',dest='lev',type=int)
parser.add_argument('-L',dest='levels',default=None,type=float,nargs=3)
parser.add_argument('-y',dest='years',default=None,help='average of this range of years. form startYear,stopYear')
parser.add_argument('-s',dest='seasons',default=['DJF','MAM','JJA','SON'],nargs='+',help='only plot those seasons')
parser.add_argument('-Y',dest='plot_years',default=None,help='Only plot these years. form startYear,stopYear')
parser.add_argument('--qbo',dest='qbo',default=None,help='Only take ensemble members which start in given QBO phase. Either + or - if given.')
args = parser.parse_args()
sns.set_context('notebook')

seasons = args.seasons
nseas = len(seasons)
olevel = args.lev
model = args.model
qmodel= fc.ModName(model)
colrs = sns.color_palette()

var = fc.variables[qmodel][args.var]

if args.plot_years is None:
    yslice = slice(None,'0010')
else:
    pyrs = ['{0:04d}'.format(int(yr)) for yr in args.plot_years.split(',')]
    yslice = slice(*pyrs)

#ctrls = xr.open_dataset(model+'_season_ctrl_ens.nc')[var].isel(time=slice(1,None))
diffs = xr.open_dataset(model+'_season_delta_ens.nc',decode_times=False)[var]#.isel(time=slice(1,None)) #
diffs,_,_ = fc.CorrectTime(diffs.to_dataset())
diffs = diffs[var].sel(time=yslice)
if args.qbo is not None:
    ds = xr.open_dataset('{0}_season_ctrl_ens.nc'.format(model),decode_times=False).isel(time=slice(1,None))
    qbo_pos,qbo_neg = fc.CheckQBO(ds,model)
    if args.qbo == '+':
        filtr = qbo_pos
    elif args.qbo == '-':
        filtr = qbo_neg
    diffs = diffs.isel(member=filtr)

dims = ac.FindCoordNames(diffs)
lev = dims['pres']
diffs = diffs.sel({lev:olevel})
#ctrls = ctrls.sel({lev:olevel})




if args.levels is None:
    levels = 20
else:
    levs = [*args.levels[:2],int(args.levels[2])]
    levs = [l for l in levs if l != 0 ]
    levels = np.linspace(*levs)

if args.years is None:
    cwrap = 2#min(5,len(diffs.time)//4)
    #if cwrap <= 5:
    #    nrows = 1
    #else:
    nrows = 5//cwrap
    filtrdim = 'time.season'
    timedim = 'time'
    ttle = ', year YY'
else:
    cwrap = 2
    nrows = int(nseas/cwrap)
    fig,axs,transf = ac.Projection('PlateCarree',ncols=cwrap,nrows=nrows,kw_args={'central_longitude':155})
    years = [int(y) for y in args.years.split(',')]
    yslce = slice('{0:04d}'.format(years[0]),'{0:04d}'.format(years[1]))
    diffs = diffs.sel(time=yslce).groupby('time.season').mean()
    filtrdim = 'season'
    timedim='season'
    ttle = ', years {0}-{1}'.format(*years)
    
    
for s,season in enumerate(seasons):
    if args.years is None:
        fig,axs,transf = ac.Projection('PlateCarree',ncols=cwrap,nrows=nrows,kw_args={'central_longitude':155})
    transf['add_colorbar'] = False
    if s == 0 or args.years is None:
        fig.set_figheight(nrows*3)
        fig.set_figwidth(cwrap*6)
    filtr = diffs[filtrdim] == season
    tmp = diffs.isel({timedim:filtr})
    pval = ac.StatTest(tmp,0,'T','member',parallel=True)
    tmp = tmp.mean('member').where(pval < 0.1)
    for a,ax in enumerate(axs.flat):
        if args.years is None:
            tmpp = tmp.isel(time=a)
            yr = tmpp.time.dt.year.values
            attle = ttle.replace('YY',str(yr))
        else:
            if s == a:
                tmpp = tmp.sel(season=season)
                attle = ttle
            else:
                continue
        cf = tmpp.plot.contourf(ax=ax,x='lon',levels=levels,extend='both',**transf)
        if not ax.get_subplotspec().is_last_row():
            ax.set_xlabel('')
        if not ax.get_subplotspec().is_first_col():
            ax.set_ylabel('')
        ax.set_title('{0}{1}, {2}{3}'.format(args.var,olevel,season,attle))
        ax.gridlines()
        ax.coastlines()
    ac.AddPanelLabels(axs,'upper left',ypos=1.12) 
    if s == len(seasons)-1 or args.years is None:
        ac.AddColorbar(fig,axs,cf)
        #fg.fig.suptitle(season)
        if args.years is None:
            suff = season
        else:
            suff = '{0}-{1}'.format(*years)
        outFile = 'figures/{3}_{0}{1}_{2}.pdf'.format(args.var,olevel,suff,model)
        if args.qbo:
            outFile = fc.RenameQBOFile(outFile,args.qbo)
        fig.savefig(outFile,bbox_inches='tight')
        print(outFile)

    
