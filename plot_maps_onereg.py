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
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',help='Choose model/run to plot.')
args = parser.parse_args()
model = args.model
sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()

fields = ['OLR','P','SLP','TS']

var_fields = {'TS':'TREFHT'}

scales = {'TS': 1, 'OLR': 1}

if 'waccm' in model:
    scales['P'] = 1e3*86400
    scales['SLP'] = 0.01
else:
    scales['P'] = 86400
    scales['SLP'] = 1

areas = {
    'TS':{
        'DJF': {
            #'Scandinavia':  {'lon':slice(10,50)  ,'lat':slice(58,70)},
            #'Eurasia': {'lon':slice(40,80)  ,'lat':slice(35,50)},
            'NAmerica':{'lon':slice(235,265),'lat':slice(45,65)},
            },
        'JJA': {
            #'Scandinavia':  {'lon':slice(20,60)  ,'lat':slice(55,70)},
            #'NAmerica': {'lon':slice(260,290),'lat':slice(40,60)},
            'Australia':{'lon':slice(120,145),'lat':slice(-28,-18)},
            },
        },
}

loncents = {
    'P' : 155,
    'OLR':155,
    'TS': 155,
    'SLP':155,
}

avg_years = {
    'P':   [2,8],
    'OLR': [2,8],
    'TS':  [2,8],
    'SLP': [2,8],
}

vmaxs = {
    'TS': 1.5,
    'P' : 1.0,
    'SLP': 6,
    'OLR': 5,
    }

cmaps = {
    'TS': 'RdBu_r',
    'P' : 'PuOr',
    'OLR':'PiYG',
    'SLP':'BrBG_r',
    }

labls = {
    'TS': 'T2m [K]',
    'P' : 'Q [mm/day]',
    'OLR':'OLR [W/m2]',
    'SLP':'SLP [hPa]',
    }

pert = xr.open_dataset('{0}_pert_ens.nc'.format(model),decode_times=False)
ctrl = xr.open_dataset('{0}_ctrl_ens.nc'.format(model),decode_times=False)


dTS = pert - ctrl
dTS,_,_ = fc.CorrectTime(dTS)

nvars = len(fields)
fig = plt.figure(figsize=[6*4,3*nvars])
transf = {'transform':ccrs.PlateCarree()}

axs = []
p = 0
for f,field in enumerate(fields):
    var = fc.variables[model][field]
    da = dTS[var]*scales[field]
    # first, seasonal maps
    for s,season in enumerate(areas[field].keys()):
        p += 1
        ax = fig.add_subplot(nvars,4,p,projection=ccrs.EqualEarth(central_longitude=loncents[field]))
        axs.append(ax)
        if isinstance(season,int):
            filtr = da['time.month'] == season
            seas = calendar.month_abbr[season]
        else:
            filtr = da['time.season'] == season
            seas = season
        dm = da.isel(time=filtr)
        dm = dm.sel(time=slice('{0:04d}'.format(avg_years[field][0]),'{0:04d}'.format(avg_years[field][1]))).mean('time')
        pval = ac.StatTest(dm,0,'T','member',parallel=True)
        dmm = dm.mean('member')
        dmm.where(pval<0.1).plot(ax=ax,vmax=vmaxs[field],cmap=cmaps[field],cbar_kwargs={'label':labls[field]},**transf)
        ax.gridlines()
        ax.coastlines()
        geoms = []
        r=0
        for region,sels in areas[field][season].items():
            if sels['lon'].start > 180:
                lonstart = sels['lon'].start - 360
            else:
                lonstart = sels['lon'].start
            if sels['lon'].stop > 180:
                lonstop = sels['lon'].stop - 360
            else:
                lonstop = sels['lon'].stop
            geom = geometry.box(minx=lonstart,maxx=lonstop,
                                miny=sels['lat'].start,maxy=sels['lat'].stop)
            #geoms.append(geom)
            ax.add_geometries([geom],edgecolor=colrs[r],facecolor='none',crs=ccrs.PlateCarree())
            r += 1
        ax.set_title(field+' anomalies, {0}, years {1}-{2}'.format(seas,avg_years[field][0],avg_years[field][1]))
    # then, regional bars
    if len(areas[field][season].items())>0:
        for season in areas[field].keys():
            p += 1
            if isinstance(season,int):
                seas = calendar.month_abbr[season]
            else:
                seas = season
            ax = fig.add_subplot(nvars,4,p)
            axs.append(ax)
            dt = []
            for region,sels in areas[field][season].items():
                tmp = ac.GlobalAvgXr(da.sel(sels)).mean('lon')
                labl = ' '.join([region,seas])
                tmp = tmp.groupby('time.year').mean()
                tmp['region'] = labl
                dt.append(tmp)
            dx = xr.concat(dt,'region')
            dx.name = var
            df = fc.MakeDataFrame(dx,'region')
            sns.barplot(data=df,x='year',y=var,hue='region',ci=90,ax=ax)
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax.grid(True,axis='x',which='minor')
            sns.despine(bottom=True,left=True)
            ax.set_title('{1} anomalies, {0}'.format(seas,field))
            ax.set_ylabel(var)
            if f == nvars-1:
                ax.set_xlabel('time [years]')
            else:
                ax.set_xlabel('')
    else:
        p += 2
ac.AddPanelLabels(axs,'upper left',ypos=1.1) 
outFile = 'figures/{0}_maps_all.pdf'.format(model)
fig.savefig(outFile,bbox_inches='tight',transparent=True)
print(outFile)
