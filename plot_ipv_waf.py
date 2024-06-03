
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
parser.add_argument('-y',dest='years',default=None,help='average of this range of years. form startYear,stopYear. Plot every year independendly of of form startYear-stopYear' )
parser.add_argument('-v',dest='variable',default='ipv',help='which variable to show in background shading.')
parser.add_argument('-l',dest='labels',default=0,type=int,help='Start labelling panels from this number (inclusive).')
parser.add_argument('--sig',dest='sig',action='store_false',help='Do NOT check for significance.')
parser.add_argument('--waf',dest='waf',action='store_false',help='Do NOT plot WAF.')
args = parser.parse_args()
model = args.model
qmodel = fc.ModName(model)
sns.set_context('notebook')
sns.set_style('whitegrid')
colrs = sns.color_palette()
var = args.variable
if var == 'ipv':
    tvar = '350K PV'
    #units= '1e6 m$^2$K/kg.s'
    units= 'PV [PVU]'
else:
    tvar = var
    units= None

seasons = args.seasons
nseas = len(seasons)

levels = {'ipv':[l for l in np.linspace(-2e-1,2e-1,21) if l !=0]}

latsel = {'lat':slice(-80,80)}
istep = 4

remove_tropics = True

dipvs = xr.open_dataset(model+'_{0}_season_delta_ens.nc'.format(var),decode_times=False)
if units is None:
    if 'units' in dipvs[var].attrs:
        units = '{0} [{1}]'.format(var,dipvs[var].attrs['units'])
    else:
        units = ''
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


time_mean = True
if args.years is not None:
    if ',' in args.years: # time mean
        yrs = [int(y) for y in args.years.split(',')]
        tslice = {'time':slice('{:04d}'.format(yrs[0]),'{:04d}'.format(yrs[1]))}
        dipvsm = dipvs.sel(tslice)
        wafsm = wafs.sel(tslice) 
    else:
        yrs = [int(y) for y in args.years.split('-')]
        time_mean = False
        dipvsm = dipvs
        wafsm = wafs

if time_mean:
    dipvsm = dipvsm.groupby('time.season').mean(keep_attrs=True)
    wafsm = wafsm.groupby('time.season').mean()
    ncols = nseas
    nrows = 1
    year_it = [0]
else:
    # shift so that first DJF is really JF of year 0
    dipvsm = dipvsm.shift(time=1).resample(time='QS-JAN').mean(keep_attrs=True)
    wafsm = wafsm.shift(time=1).resample(time='QS-JAN').mean()
    ncols = 2
    nyears = yrs[1]-yrs[0]+1
    year_it = range(yrs[0],yrs[1]+1)
    nrows = int(np.ceil(nyears/ncols))

fig,axs,transf = ac.Projection('PlateCarree',ncols=ncols,nrows=nrows,kw_args={'central_longitude':155})
fig.set_figheight(3*nrows)
fig.set_figwidth(6*ncols)

if args.sig:
    pval = ac.StatTest(dipvsm,0,'T','member',parallel=True)
    dipvsm = dipvsm.mean('member',keep_attrs=True).where(pval<0.1)
else:
    dipvsm = dipvsm.mean('member',keep_attrs=True).sel(latsel)
dipvsm = ac.CloseGlobe(dipvsm)

wafa = np.sqrt(wafsm.wx**2+wafsm.wy**2)
if args.sig:
    pval = ac.StatTest(wafa,0,'T','member',parallel=True)
    wafsm = wafsm.mean('member').where(pval<0.1)
else:
    wafsm = wafsm.mean('member')

p = -1
for y,year in enumerate(year_it):
    for s,seas in enumerate(seasons):
        p += 1
        ax = axs.flatten()[p]
        if 'season' in dipvsm.dims:
            filtrd = dipvsm.sel(season=seas)
            filtrw = wafsm.sel(season=seas)
        else:
            filtr = (dipvsm.time.dt.season == seas)*(dipvsm.time.dt.year == year)
            filtrd= dipvsm.isel(time=filtr).squeeze()
            filtr = (wafsm.time.dt.season == seas)*(wafsm.time.dt.year == year)
            filtrw = wafsm.isel(time=filtr).squeeze()
        cf=filtrd.plot.contourf(levels=levels[var],ax=ax,add_colorbar=False,extend='both',**transf)
        if args.waf:
            Q = filtrw.plot.quiver(x='lon',y='lat',u='wx',v='wy',ax=ax,scale=3e2,**transf)
            if p == 0:
                ax.quiverkey(Q,1.1,0.5,20,r'20 $\frac{\mathrm{m^2}}{\mathrm{s^2}}$')
        ax.gridlines()
        ax.coastlines()
        ttle = tvar
        if args.waf:
            ttle = ttle + ' & 200hPa WAF'
        ttle = ttle + ', {0}'.format(seas)
        if not time_mean:
            ttle = ttle.replace(seas,'year {0} {1}'.format(year,seas))
        ax.set_title(ttle)
        ax.set_extent([0,360,-90,90],crs=transf['transform'])
ac.AddColorbar(fig,axs,cf,cbar_args={'extend':'both','label':units})
ac.AddPanelLabels(axs,'upper left',ypos=1.15,start_index=args.labels)
if isinstance(args.seasons,list):
    if len(args.seasons)>1:
        seasarg = '_'.join(args.seasons)
    else:
        seasarg = args.seasons[0]
else:
    seasarg = args.seasons
outFile = 'figures/{0}_{1}_waf_{2}_year.pdf'.format(model,var,seasarg)
if not args.waf:
    outFile = outFile.replace('_waf','')
if args.years is not None:
    outFile = outFile.replace('_year.pdf','_years{0}-{1}.pdf'.format(*yrs))
if args.qbo is not None:
    outFile = fc.RenameQBOFile(outFile,args.qbo)
fig.savefig(outFile,bbox_inches='tight',transparent=True)
print(outFile)

