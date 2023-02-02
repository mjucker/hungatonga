#!env python
import xarray as xr
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
from matplotlib import pyplot as plt
import calendar
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',help='Choose model/run to plot.')
args = parser.parse_args()
model = args.model
sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()

max_year = 10
avg_years= [2,10]
vmax = 6.0e0
scale = 0.01

var = 'SLP'

pert = xr.open_dataset('{0}_pert_ens.nc'.format(model),decode_times=False)#.TS
ctrl = xr.open_dataset('{0}_ctrl_ens.nc'.format(model),decode_times=False)#.TS

pert = pert[fc.variables[model][var]]
ctrl = ctrl[fc.variables[model][var]]

dTS = pert - ctrl
dTS,_,_ = fc.CorrectTime(dTS.to_dataset())
dTS = dTS[fc.variables[model][var]]
dTS = dTS.sel(time=slice(None,'{0:04d}'.format(max_year)))
dTS = scale*dTS

#dTS2 = fc.ReadFile('waccm_pert_ens.nc').TS - fc.ReadFile('waccm_ctrl_ens.nc').TS


#select_seasons = ['DJF','JJA']
select_months = [1,7]

#months = {'Eurasia': 'DJF',
#          'NAmerica':'DJF',
#          #'Arctic' : 'DJF',
#          #'EAsia':'DJF',
#          #'Australia': 'DJF',
#}
areas = {7:{
                'ASL':  {'lon':slice(235,275)  ,'lat':slice(-70,-60)},
               },
         1:{
                'AL':  {'lon':slice(180.1,220)  ,'lat':slice(40,60)},
               },
         #'Arctic' : {'lon':slice(0,360),  'lat':slice(80,90)},
         #'EAsia':   {'lon':slice(100,125),'lat':slice(40,55)},
         #'Australia':{'lon':slice(115,155),'lat':slice(-38,-20)}
         }

from shapely import geometry
from cartopy import crs as ccrs

for select_month in select_months:
    mam = dTS['time.month'] == select_month
    dTS_MAM = dTS.isel(time=mam).sel(time=slice('{0:04d}'.format(avg_years[0]),'{0:04d}'.format(avg_years[1]))).mean('time')
    fig,ax,transf = ac.Projection('PlateCarree',kw_args={'central_longitude':180})
    fig.set_figwidth(6)
    pval = ac.StatTest(dTS_MAM,0,'T','member',parallel=True)
    dTS_MAM.mean('member').where(pval<0.1).plot(ax=ax,vmax=vmax,cbar_kwargs={'shrink':0.5},**transf)
    ax.gridlines()
    ax.coastlines()
    geoms = []
    for region,sels in areas[select_month].items():
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
        geoms.append(geom)
    ax.add_geometries(geoms,edgecolor='k',facecolor='none',crs=ccrs.PlateCarree())
    ax.set_title(var+' anomalies, {0}, years {1}-{2}'.format(calendar.month_abbr[select_month],avg_years[0],avg_years[1]))
    outFile = 'figures/{1}_{2}_{0}.pdf'.format(calendar.month_abbr[select_month],model,var)
    fig.savefig(outFile,transparent=True,bbox_inches='tight')
    print(outFile)
    #
    #
    dt = []
    for region,sels in areas[select_month].items():
        tmp = ac.GlobalAvgXr(dTS.sel(sels)).mean('lon')
        filtr = tmp['time.month'] == select_month
        labl = region+' {0}'.format(calendar.month_abbr[select_month])
        #mnx= dTS.sel(sels).mean('member').isel(time=filtr).groupby('time.year').mean().sel(year=slice(2,None)).mean('year')
        mnx= dTS.sel(sels).mean('member').sel(time=slice('{0:04d}-12'.format(avg_years[0]),'{0:04d}'.format(avg_years[1]))).groupby('time.month').mean()
        mn = mnx.sel(month=select_month).min().values
        mx = mnx.sel(month=select_month).max().values
        print('{0}: min = {1:+.2f}, max = {2:+.2f}'.format(labl,mn,mx))
        tmp = tmp.isel(time=filtr).groupby('time.year').mean()
        tmp['region'] = labl
        dt.append(tmp)
    dx = xr.concat(dt,'region')
    dx.name = 'Ts'
    df = fc.MakeDataFrame(dx,'region')

    fig,ax = plt.subplots(figsize=[6,3])
    sns.barplot(data=df,x='year',y='Ts',hue='region',ci=90,ax=ax)
    sns.despine(bottom=True,left=True)
    ax.set_title('{1} anomalies, {0}'.format(calendar.month_abbr[select_month],var))
    ax.set_ylabel(var)
    ax.set_xlabel('time [years]')
    outFile = 'figures/{0}_{2}_regions_{1}.pdf'.format(model,calendar.month_abbr[select_month],var)
    fig.savefig(outFile,bbox_inches='tight',transparent=True)
    print(outFile)
