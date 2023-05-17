#!env python
import xarray as xr
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from shapely import geometry
from cartopy import crs as ccrs
import calendar
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',help='Choose model/run to plot.')
parser.add_argument('-v',dest='vars',default=None,nargs='+')
args = parser.parse_args()
model = args.model
sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()

fields = ['OLR','P','SLP','TS']
if args.vars is not None:
    fields = args.vars

rolls = [1,3,6,12]#,24,36,48,60]
yrolls= [1]#,2,3,4,5]

var_fields = {'TS':'TREFHT'}

areas = {
    'TS':{
        
        'Scandinavia':  {'lon':slice(10,50)  ,'lat':slice(58,70)},
        'Eurasia': {'lon':slice(40,80)  ,'lat':slice(35,50)},
        'NAmerica':{'lon':slice(235,265),'lat':slice(45,65)},
        'Australia':{'lon':slice(120,145),'lat':slice(-28,-18)},
    },
    'SLP':{
        'ASL':  {'lon':slice(235,275)  ,'lat':slice(-70,-60)},
        'AL':  {'lon':slice(180.1,220)  ,'lat':slice(40,60)},
    },
    'P': {
        'ITCZ':{'lon':slice(180.1,210), 'lat': slice(5,10)},
        'MC'  :{'lon':slice(120,140),   'lat': slice(0,20)},
        'IO':  {'lon':slice(40,95),     'lat': slice(-12,0)},
        'Australia':{'lon':slice(120,145),'lat':slice(-28,-18)},
    },
    'OLR':{
        'ITCZ':{'lon':slice(180.1,210), 'lat': slice(5,10)},
        'MC'  :{'lon':slice(120,140),   'lat': slice(0,20)},
        'IO':  {'lon':slice(40,95),     'lat': slice(-12,0)},
        'Australia':{'lon':slice(120,145),'lat':slice(-28,-18)},
    },
}


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

# plot a map with all boxes:
unique_regions = {}
for var in areas.keys():
    for region in areas[var].keys():
        if region not in unique_regions.keys():
            unique_regions[region] = areas[var][region]
fig,ax,transf = ac.Projection('Robinson',kw_args={'central_longitude':155})
for region,sels in unique_regions.items():
    fc.AddBox(ax,sels,'r')
ax.coastlines()
outFile = 'figures/map_of_regions.pdf'
fig.savefig(outFile,transparent=True,bbox_inches='tight')
print(outFile)


# now start working with the actual data
pert = xr.open_dataset('{0}_pert_ens.nc'.format(model),decode_times=False)
ctrl = xr.open_dataset('{0}_ctrl_ens.nc'.format(model),decode_times=False)

CT,_,_ = fc.CorrectTime(ctrl)

dTS = pert - ctrl
dTS,_,_ = fc.CorrectTime(dTS)

# analyze monthly timeseries
for s,roll in enumerate(rolls):
    for f,field in enumerate(areas.keys()):
        axs = []
        nreg = len(areas[field].keys())
        ncols = 2
        nrows = (nreg)//ncols
        fig = plt.figure(figsize=[ncols*4,nrows*3])
        for r,reg in enumerate(areas[field].keys()):
            ax = fig.add_subplot(nrows,ncols,r+1)
            axs.append(ax)
            var = fc.variables[model][field]
            if var in var_fields.keys():
                var = var_fields[var]
            da = dTS[var]*scales[field]
            dc = CT[var]*scales[field]
            tmp = da.sel(areas[field][reg]).mean(['lon','lat'])
            tmp = tmp.rolling(time=roll).mean()
            std = dc.sel(areas[field][reg]).mean(['lon','lat']).rolling(time=roll).mean().std('member')
            tmp = tmp/std
            tmp.plot.line(x='time',color=colrs[r],alpha=0.1,add_legend=False,ax=ax)
            tmp.mean('member').plot.line(ax=ax,color=colrs[r],lw=1,ls='--',add_legend=False)
            pval = ac.StatTest(tmp,0,'T','member')
            tmp.mean('member').where(pval<0.1).plot.line(ax=ax,color=colrs[r],lw=2,add_legend=False)
            ax.set_title(reg)
            if r < 2:
                ax.set_xlabel('')
        fig.suptitle('{0}, {1}-month rolling mean'.format(labls[field],roll))
        ac.AddPanelLabels(axs,'upper left',ypos=1.1) 
        outFile = 'figures/{0}_{1}_roll{2}_lines_regions.pdf'.format(model,field,roll)
        fig.savefig(outFile,bbox_inches='tight',transparent=True)
        print(outFile)

# now do the same thing by season
sTS = xr.open_dataset('waccm_season_delta_ens.nc')
sCT = xr.open_dataset('waccm_season_ctrl_ens.nc')
#sTS = dTS.resample(time='QS-DEC').mean()
#sCT = CT.resample(time='QS-DEC').mean()

seasons = np.unique(sTS.time.dt.season)
# plot per variable per roll per region, 1 panel for each season
for f,field in enumerate(areas.keys()):
    var = fc.variables[model][field]
    if var in var_fields.keys():
        var = var_fields[var]
    da = sTS[var]*scales[field]
    dc = sCT[var]*scales[field]
    for a,roll in enumerate(yrolls):
        for r,reg in enumerate(areas[field].keys()):
            axs = []
            ncols = 2
            nrows = 2
            fig = plt.figure(figsize=[ncols*4,nrows*3])
            for s,seas in enumerate(seasons):
                ax = fig.add_subplot(nrows,ncols,s+1)
                axs.append(ax)
                filtr = sTS[var].time.dt.season == seas
                tmp = da.sel(areas[field][reg]).mean(['lon','lat']).isel(time=filtr)
                tmp = tmp.rolling(time=roll).mean()
                std = dc.sel(areas[field][reg]).mean(['lon','lat']).isel(time=filtr)
                std = std.rolling(time=roll).mean().std('member')
                tmp = tmp/std
                tmp.plot.line(x='time',color=colrs[r],alpha=0.1,add_legend=False,ax=ax)
                tmp.mean('member').plot.line(ax=ax,color=colrs[r],lw=1,ls='--',add_legend=False)
                pval = ac.StatTest(tmp,0,'T','member')
                tmp.mean('member').where(pval<0.1).plot.line(ax=ax,color=colrs[r],lw=2,add_legend=False)
                ax.set_title(seas)
                if s < 2:
                    ax.set_xlabel('')
            fig.suptitle('{0}, {1}, {2}-year rolling mean'.format(labls[field],reg,roll))
            ac.AddPanelLabels(axs,'upper left',ypos=1.1) 
            outFile = 'figures/{0}_{1}_{2}_roll{3}_lines.pdf'.format(model,field,reg,roll)
            fig.savefig(outFile,bbox_inches='tight',transparent=True)
            print(outFile)
