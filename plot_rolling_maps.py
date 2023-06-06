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
parser.add_argument('--qbo',dest='qbo',default=None,help='Only take ensemble members which start in given QBO phase. Either + or - if given.')
args = parser.parse_args()
model = args.model
qmodel = fc.ModName(model)
sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()

if args.vars is None:
    fields = ['OLR','P','SLP','TS']
    if model.lower() == 'waccm':
        fields.append('CLDTOT')
else:
    fields = args.vars

yrolls = [1]
#yrolls = [2,3]

var_fields = {'waccm':{'TS':'TREFHT'},
              'mima' : {},
              'hthh':{},
              'hthh_fix':{},
              'bench_SH':{},
              }

cmap_vars = {'TS':'t','P':'precip','OLR':'olr','SLP':'slp','default':'RdBu_r','CLDTOT': 'cldfrac'}


scales = {'TS': 1, 'OLR': 1, 'CLDTOT': 1}

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
    'CLDTOT':'CLDTOT [.]',
    }

# note: it is not practical to plot maps on the monthly rolling means
#        as this would create one figure per field per rolling mean with one panel per month
#        in the simulation, which is up to 120 months!

# create maps by season
sTS = xr.open_dataset('{0}_season_delta_ens.nc'.format(model),decode_times=False).isel(time=slice(1,None))
sCT = xr.open_dataset('{0}_season_ctrl_ens.nc'.format(model),decode_times=False).isel(time=slice(1,None))
sTS,_,_ = fc.CorrectTime(sTS)
sCT,_,_ = fc.CorrectTime(sCT)

if args.qbo is not None:
    U = fc.variables[model]['U']
    lev = ac.FindCoordNames(sCT)['pres']
    U50 = sCT[U].sel({lev:50,'lat':slice(-5,5)}).isel(time=0).mean(['lon','lat'])
    del U50['time']
    del U50[lev]
    qbo_pos = U50 > 0
    qbo_neg = U50 < 0

#sTS = dTS.resample(time='QS-DEC').mean()
#sCT = CT.resample(time='QS-DEC').mean()

seasons = np.unique(sTS.time.dt.season)
# plot per variable per roll per season, 1 panel for each rolling year
for f,field in enumerate(fields):
    var = fc.variables[qmodel][field]
    if var in var_fields[qmodel].keys():
        var = var_fields[qmodel][var]
    if field in cmap_vars.keys():
        cmap = at.cmaps[cmap_vars[field]]
    else:
        cmap = cmap_vars['default']
    da = sTS[var]*scales[field]
    dc = sCT[var]*scales[field]
    if args.qbo is not None:
        if args.qbo == '+':
            filtr = qbo_pos
        elif args.qbo == '-':
            filtr = qbo_neg
        da = da.isel(member=filtr)
        dc = dc.isel(member=filtr)
    for a,roll in enumerate(yrolls):
        for s,seas in enumerate(seasons):
            filtr = sTS[var].time.dt.season == seas
            if roll == 1:
                tmp = da.isel(time=filtr)
                std = dc.isel(time=filtr).std('member')
            else:
                tmp = da.isel(time=filtr).rolling(time=roll).mean()
                std = dc.isel(time=filtr).rolling(time=roll).mean().std('member')
            ntimes = len(tmp.time.values)
            ncols = min(ntimes,5)
            nrows = (ntimes-1)//ncols+1
            fig,axs,transf = ac.Projection('PlateCarree',ncols=ncols,nrows=nrows,kw_args={'central_longitude':155})
            fig.set_figheight(nrows*3)
            fig.set_figwidth(ncols*4)
            tmp = tmp/std
            #pval = ac.StatTest(tmp,0,'T','member',parallel=True)
            #filtr = pval < 0.1
            pval = ac.StatTest(tmp,tmp.mean('member'),'sign','member',parallel=True)
            filtr = pval > 0.66 # at least 20 members have same sign as mean
            for t,time in enumerate(tmp.time):
                ax = axs.flatten()[t]
                tmp.mean('member').where(filtr).isel(time=t).plot.contourf(levels=20,ax=ax,cmap=cmap,robust=True,**transf)
                ax.coastlines()
                ax.set_title('year {0}'.format(time.dt.year.values))
            ttle = '{0}, {1}'.format(labls[field],seas)
            outFile= 'figures/{0}_{1}_{2}_maps.pdf'.format(model,field,seas)
            if roll > 1:
                ttle = ttle + ', {0}-year rolling mean'.format(roll)
                outFile=outFile.replace('_maps.pdf','_roll{0}_maps.pdf'.format(roll))
            if args.qbo is not None:
                ttle = ttle+', QBO'+args.qbo
                qnme = {'+':'p','-':'m'}
                outFile=outFile.replace('.pdf','_QBO{0}.pdf'.format(qnme[args.qbo]))
            fig.suptitle(ttle)
            ac.AddPanelLabels(axs,'upper left',ypos=1.2) 
            fig.savefig(outFile,bbox_inches='tight',transparent=True)
            print(outFile)
        plt.close('all')
