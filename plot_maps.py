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

trans_fields = {'TS':'TREFHT'}

scales = {'TS': 1, 'OLR': 1}

if 'waccm' in model:
    scales['P'] = 1e3*86400
    scales['SLP'] = 0.01
else:
    scales['P'] = 86400
    scales['SLP'] = 1

seasons = {
    'TS': ['DJF','JJA'],
    'OLR':['DJF','JJA'],
    'P'  :['DJF','JJA'],
    'SLP':[1,7],
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
    'OLR': 5.0,
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

pert = xr.open_dataset('{0}_pert_ens.nc'.format(model),decode_times=False)
ctrl = xr.open_dataset('{0}_ctrl_ens.nc'.format(model),decode_times=False)


dTS = pert - ctrl
dTS,_,_ = fc.CorrectTime(dTS)

nvars = len(fields)
fig = plt.figure(figsize=[6*2,3*nvars])
transf = {'transform':ccrs.PlateCarree()}

axs = []
p = 0
for f,field in enumerate(fields):
    var = fc.variables[model][field]
    if var in trans_fields.keys():
        var = trans_fields[var]
    da = dTS[var]*scales[field]
    # first, seasonal maps
    for s,season in enumerate(seasons[field]):
        p += 1
        ax = fig.add_subplot(nvars,2,p,projection=ccrs.Robinson(central_longitude=loncents[field]))
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
        dmm.where(pval<0.1).plot(ax=ax,vmax=vmaxs[field],cmap=cmaps[field],cbar_kwargs={'label':labls[field],'shrink':0.85},**transf)
        ax.gridlines()
        ax.coastlines()
        if field in areas.keys():
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
                ax.add_geometries([geom],edgecolor=colrs[r],facecolor='none',crs=ccrs.PlateCarree(),zorder=5)
                r += 1
        ax.set_title(field+' anomalies, {0}, years {1}-{2}'.format(seas,avg_years[field][0],avg_years[field][1]))
ac.AddPanelLabels(axs,'upper left',ypos=1.1) 
outFile = 'figures/{0}_maps_boxed.pdf'.format(model)
fig.savefig(outFile,bbox_inches='tight',transparent=True)
print(outFile)

# now add bar charts for the regions
for field in areas.keys():
    var = fc.variables[model][field]
    if var in trans_fields.keys():
        var = trans_fields[var]
    da = dTS[var]*scales[field]
    fig,ax = plt.subplots()
    dt = []
    for season in areas[field].keys():
        filtr = da['time.season'] == season
        das = da.isel(time=filtr)
        for region,sels in areas[field][season].items():
            tmp = ac.GlobalAvgXr(das.sel(sels)).mean('lon')
            labl = ' '.join([region,season])
            tmp = tmp.groupby('time.year').mean()
            tmp['region'] = labl
            dt.append(tmp)
    dx = xr.concat(dt,'region')
    dx.name = field
    df = fc.MakeDataFrame(dx,'region')
    sns.barplot(data=df,x='year',y=field,hue='region',ci=90,ax=ax)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(True,axis='x',which='minor')
    sns.despine(bottom=True,left=True)
    ax.set_title('{0} anomalies'.format(field))
    ax.set_ylabel(labls[field])
    ax.set_xlabel('year')
    outFile = 'figures/{0}_bars_{1}.pdf'.format(model,field)
    fig.savefig(outFile,bbox_inches='tight',transparent=True)
    print(outFile)
    
            
