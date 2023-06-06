#!env python
import xarray as xr
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
from matplotlib import pyplot as plt
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',help='Choose model/run to plot.')
parser.add_argument('-s',dest='savefig',action="store_false",help='do NOT save the figures.')
parser.add_argument('--qbo',dest='qbo',default=None,help='Only take ensemble members which start in given QBO phase. Either + or - if given.')
args = parser.parse_args()
model = args.model
sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()

max_year = 10
avg_years= [4,7]
vmax = 1.5
#vmax = None

var = 'TS'# I am interested in land temps only here
#var = 'P'
#'TREFHT'#'TS'

pert = xr.open_dataset('{0}_pert_ens.nc'.format(model),decode_times=False)#.TS
ctrl = xr.open_dataset('{0}_ctrl_ens.nc'.format(model),decode_times=False)#.TS

if args.qbo is not None:
    qbo_pos,qbo_neg = fc.CheckQBO(pert,model)
    if args.qbo == '+':
        pert = pert.isel(member=qbo_pos)
        ctrl = ctrl.isel(member=qbo_pos)
    elif args.qbo == '-':
        pert = pert.isel(member=qbo_neg)
        ctrl = ctrl.isel(member=qbo_neg)

pert = pert[fc.variables[model][var]]
ctrl = ctrl[fc.variables[model][var]]

dTS = pert - ctrl
dTS,_,_ = fc.CorrectTime(dTS.to_dataset())
dTS = dTS[fc.variables[model][var]]
dTS = dTS.sel(time=slice(None,'{0:04d}'.format(max_year)))

#dTS2 = fc.ReadFile('waccm_pert_ens.nc').TS - fc.ReadFile('waccm_ctrl_ens.nc').TS
colrs = ['r','b']
fig,ax = plt.subplots()
for r,reg in enumerate([[-60,-30],[30,60]]):
    dTm = ac.GlobalAvgXr(dTS.where(dTS!=0).mean('lon'),reg)
    pval = ac.StatTest(dTm,0,'T','member',parallel=True)
    dTmm = dTm.mean('member')
    dTms = dTm.std('member')
    if reg[0] > 0:
        hem='N'
    else:
        hem='S'
    labl = '{0}-{1}{2}'.format(abs(reg[0]),abs(reg[1]),hem)
    p = dTmm.plot(ax=ax,color=colrs[r],ls='--',lw=1)
    dTmm.where(pval<0.1).plot(ax=ax,color=colrs[r],label=labl)
    x = p[0].get_xdata()
    ax.fill_between(x,dTmm-dTms,dTmm+dTms,color=colrs[r],alpha=0.3)
    ax.legend()
#xlims = ax.get_xlim() # doesn't work with cftime
#ax.set_xlim(1,xlims[1])
#ax.set_xticks(np.arange(1,xlims[1]+1))
sns.despine(left=True,bottom=True)
ax.set_xlabel('time [years]')
ax.set_title('Land Ts anomalies [K]')
outFile = 'figures/{0}_{1}_land_hemi.pdf'.format(model,var)
if args.qbo is not None:
    outFile = fc.RenameQBOFile(outFile,args.qbo)
if args.savefig:
    fig.savefig(outFile,bbox_inches='tight',robust=True)
    print(outFile)


select_seasons = ['DJF','JJA','MAM','SON']

#months = {'Eurasia': 'DJF',
#          'NAmerica':'DJF',
#          #'Arctic' : 'DJF',
#          #'EAsia':'DJF',
#          #'Australia': 'DJF',
#}
areas = {'DJF':{
                #'Scandinavia':  {'lon':slice(10,50)  ,'lat':slice(58,70)},
                'Eurasia': {'lon':slice(40,80)  ,'lat':slice(35,50)},
                'NAmerica':{'lon':slice(235,265),'lat':slice(45,65)},
                #'Australia':{'lon':slice(120,145),'lat':slice(-28,-18)},
                'NAfrica': {'lon':slice(-10,30)  ,'lat':slice(20,35)},
               },
         'JJA':{
                #'Scandinavia':  {'lon':slice(20,60)  ,'lat':slice(55,70)},
                #'NAmerica': {'lon':slice(260,290),'lat':slice(40,60)},
                'Australia':{'lon':slice(120,145),'lat':slice(-28,-18)},
               },
         'MAM':{
                'Siberia': {'lon':slice(50,90),'lat':slice(45,65)},
               },
         'SON': {
                'NAsia'  : {'lon':slice(70,120),'lat':slice(45,55)},
               },
         #'Arctic' : {'lon':slice(0,360),  'lat':slice(80,90)},
         #'EAsia':   {'lon':slice(100,125),'lat':slice(40,55)},
         #'Australia':{'lon':slice(115,155),'lat':slice(-38,-20)}
         }

from shapely import geometry
from cartopy import crs as ccrs

for select_season in select_seasons:
    mam = dTS['time.season'] == select_season
    if select_season == 'DJF':
        dTS_MAM = dTS.isel(time=mam).sel(time=slice('{0:04d}-12'.format(avg_years[0]-1),'{0:04d}'.format(avg_years[1]))).mean('time')
    else:
        dTS_MAM = dTS.isel(time=mam).sel(time=slice('{0:04d}'.format(avg_years[0]),'{0:04d}'.format(avg_years[1]))).mean('time')
    fig,ax,transf = ac.Projection('Robinson',kw_args={'central_longitude':155})
    fig.set_figwidth(6)
    pval = ac.StatTest(dTS_MAM,0,'T','member',parallel=True)
    dTS_MAM.mean('member').where(pval<0.1).plot(ax=ax,vmax=vmax,cbar_kwargs={'shrink':0.5},**transf)
    ax.gridlines()
    ax.coastlines()
    geoms = []
    for region,sels in areas[select_season].items():
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
    ax.set_title('Ts anomalies, {0}, years {1}-{2}'.format(select_season,avg_years[0],avg_years[1]))
    outFile = 'figures/{1}_{2}_{0}.pdf'.format(select_season,model,var)
    if args.qbo is not None:
        outFile = fc.RenameQBOFile(outFile,args.qbo)
    if args.savefig:
        fig.savefig(outFile,transparent=True,bbox_inches='tight')
        print(outFile)
    #
    #
    dt = []
    for region,sels in areas[select_season].items():
        tmp = ac.GlobalAvgXr(dTS.sel(sels)).mean('lon')
        if isinstance(select_season,int):
            filtr = tmp['time.month'] == select_season
            labl = region+' month {0}'.format(select_season)
        else:
            #if select_season == 'DJF':
            #    # for DJF, we want D to be together with JF to the year after
            #    #  for this, shift data by 3 months and pretend to do MAM mean
            #    tmp = tmp.shift(time=3)
            #    filtr = tmp['time.season'] == 'MAM'
            #else:
            #    filtr = tmp['time.season'] == select_season
            filtr = tmp['time.season'] == select_season
            labl = region+' '+select_season
        #mnx= dTS.sel(sels).mean('member').isel(time=filtr).groupby('time.year').mean().sel(year=slice(2,None)).mean('year')
        mnx= dTS.sel(sels).mean('member').sel(time=slice('{0:04d}-12'.format(avg_years[0]-1),'{0:04d}'.format(avg_years[1]))).groupby('time.season').mean()
        mn = mnx.sel(season=select_season).min().values
        mx = mnx.sel(season=select_season).max().values
        print('{0}: min = {1:+.2f}K, max = {2:+.2f}K'.format(labl,mn,mx))
        tmp = tmp.isel(time=filtr).groupby('time.year').mean()
        tmp['region'] = labl
        dt.append(tmp)
    if len(dt) > 0:
        dx = xr.concat(dt,'region')
        dx.name = 'Ts'
        df = fc.MakeDataFrame(dx,'region')

        fig,ax = plt.subplots(figsize=[6,3])
        sns.barplot(data=df,x='year',y='Ts',hue='region',errorbar=('ci',90),ax=ax)
        sns.despine(bottom=True,left=True)
        ax.set_title('Ts anomalies [K], {0}'.format(select_season))
        ax.set_ylabel('Ts [K]')
        ax.set_xlabel('time [years]')
        outFile = 'figures/{0}_{2}_regions_{1}.pdf'.format(model,select_season,var)
        if args.qbo is not None:
            outFile = fc.RenameQBOFile(outFile,args.qbo)
        if args.savefig:
            fig.savefig(outFile,bbox_inches='tight',transparent=True)
            print(outFile)
