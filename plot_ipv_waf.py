import xarray as xr
import numpy as np
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
parser.add_argument('--qbo',dest='qbo',default=None,help='Only take ensemble members which start in given QBO phase. Either + or - if given.')
parser.add_argument('-s',dest='seasons',default=['DJF','MAM','JJA','SON'],nargs='+',help='only plot those seasons')
parser.add_argument('-y',dest='years',default=None,help='average of this range of years. form startYear,stopYear')
parser.add_argument('-v',dest='variable',default='ipv',help='which variable to show in background shading.')
args = parser.parse_args()
model = args.model
qmodel = fc.ModName(model)
sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()
var = args.variable
if var == 'ipv':
    tvar = '350K PV'
else:
    tvar = var

seasons = args.seasons
nseas = len(seasons)

levels = {'ipv':[l for l in np.linspace(-2e-1,2e-1,21) if l !=0]}

latsel = {'lat':slice(-80,80)}
istep = 4

remove_tropics = True

dipvs = xr.open_dataset(model+'_{0}_season_delta_ens.nc'.format(var),decode_times=False)
dipvs,_,_ = fc.CorrectTime(dipvs)
dipvs = dipvs[var]*1e6#.sel(latsel)*1e6

wafs = xr.open_dataset(model+'_waf_season_delta_ens.nc',decode_times=False).sel(lev=200)
wafs,_,_ = fc.CorrectTime(wafs)
if remove_tropics:
    wafs = wafs.where(np.abs(wafs.lat)>15)
wafs = wafs.sel(latsel).isel(lon=slice(None,None,istep),lat=slice(None,None,istep))

if args.qbo is not None:
    ds = xr.open_dataset('{0}_season_ctrl_ens.nc'.format(model),decode_times=False).isel(time=slice(1,None))
    qbo_pos,qbo_neg = fc.CheckQBO(ds,qmodel)
    if args.qbo == '+':
        filtr = qbo_pos
    elif args.qbo == '-':
        filtr = qbo_neg
    dipvs = dipvs.isel(member=filtr)
    wafs = wafs.isel(member=filtr)


if args.years is not None:
    yrs = [int(y) for y in args.years.split(',')]
    tslice = {'time':slice('{:04d}'.format(yrs[0]),'{:04d}'.format(yrs[1]))}
    dipvsm = dipvs.sel(tslice).groupby('time.season').mean()
    wafsm = wafs.sel(tslice).groupby('time.season').mean()

    fig,axs,transf = ac.Projection('PlateCarree',ncols=nseas,kw_args={'central_longitude':155})
    fig.set_figheight(3)
    fig.set_figwidth(6*nseas)

    pval = ac.StatTest(dipvsm,0,'T','member',parallel=True)
    dipvsm = dipvsm.mean('member').where(pval<0.1)
    dipvsm = ac.CloseGlobe(dipvsm)

    wafa = np.sqrt(wafsm.wx**2+wafsm.wy**2)
    pval = ac.StatTest(wafa,0,'T','member',parallel=True)
    wafsm = wafsm.mean('member').where(pval<0.1)

    for s,seas in enumerate(seasons):
        cf=dipvsm.sel(season=seas).plot.contourf(levels=levels[var],ax=axs[s],add_colorbar=False,extend='both',**transf)
        wafsm.sel(season=seas).plot.quiver(x='lon',y='lat',u='wx',v='wy',ax=axs[s],**transf)
        axs[s].gridlines()
        axs[s].coastlines()
        axs[s].set_title(tvar+' & 200hPa WAF, {0}'.format(seas))
        axs[s].set_extent([0,360,-90,90],crs=transf['transform'])
    ac.AddColorbar(fig,axs,cf,cbar_args={'extend':'both'})
    ac.AddPanelLabels(axs,'upper left',ypos=1.15)
    outFile = 'figures/{0}_{1}_waf_season_year.pdf'.format(model,var)
    if args.years is not None:
        outFile = outFile.replace('_year.pdf','_years{0}-{1}.pdf'.format(*yrs))
    if args.qbo is not None:
        outFile = fc.RenameQBOFile(outFile,args.qbo)
    fig.savefig(outFile,bbox_inches='tight',transparent=True)
    print(outFile)

