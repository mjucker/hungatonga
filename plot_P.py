#!env python
import xarray as xr
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
from matplotlib import pyplot as plt
import calendar
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',help='Choose model/run to plot.')
args = parser.parse_args()
model = args.model

sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()

variables = fc.variables

PREC = variables[model]['P']
OLR  = variables[model]['OLR']

zoom = {'lon':slice(50,280),'lat':slice(-60,60)}
zoom = None

avg_years = [1,5]

vmaxs = {'PREC':1.0, 'FLNT': 5.0}

if model == 'waccm':
    to_mm_per_day = 1e3*86400
else:
    to_mm_per_day = 86400

tmp = xr.open_dataset('{0}_pert_ens.nc'.format(model),decode_times=False)
pert = xr.merge([tmp[PREC]*to_mm_per_day,tmp[OLR]])
tmp = xr.open_dataset('{0}_ctrl_ens.nc'.format(model),decode_times=False)
ctrl = xr.merge([tmp[PREC]*to_mm_per_day,tmp[OLR]])

dTS = pert - ctrl
dTS,_,_ = fc.CorrectTime(dTS)
dTS = dTS.sel(zoom)

#dTS2 = fc.ReadFile('waccm_pert_ens.nc').TS - fc.ReadFile('waccm_ctrl_ens.nc').TS
select_month = 1
select_seasons = ['DJF','JJA']
max_year = 10

dTS = dTS.sel(time=slice(None,'{0:04d}'.format(max_year)))

cmaps = {OLR:'PiYG',PREC:'PuOr'}
if OLR != 'OLR':
    cmaps['OLR'] = cmaps[OLR]
cmaps['GPCP']= cmaps[PREC]

mam = dTS['time.month'] == select_month

hatch = ['///']

labels = {'FLNT':'OLR [W/m2]','PREC':'Q [mm/day]'}

mins = {OLR:-15,PREC:-2.5}
if OLR != 'OLR':
    mins['OLR'] = 3*mins[OLR]
mins['GPCP']= 3*mins[PREC]

for var in dTS.data_vars:
    dTS_MAM = dTS[var].isel(time=mam).sel(time=slice('{0:04d}'.format(avg_years[0]),'{0:04d}'.format(avg_years[1]))).mean('time')
    fig,ax,transf = ac.Projection('PlateCarree',kw_args={'central_longitude':180})
    transf['cbar_kwargs'] = {'shrink':0.9,'label':labels[var]}
    fig.set_figwidth(12)
    pval = ac.StatTest(dTS_MAM,0,'T','member',parallel=True)
    dTS_MAM.mean('member').plot(ax=ax,cmap=cmaps[var],**transf)
    del transf['cbar_kwargs']
    dTS_MAM.mean('member').where(pval<0.1).plot.contourf(ax=ax,colors='none',hatches=hatch,add_colorbar=False,**transf)
    ax.gridlines()
    ax.coastlines()
    ax.set_title('{0} anomalies, {1}, years {2}-{3}'.format(var,calendar.month_abbr[select_month],avg_years[0],avg_years[1]))
    outFile = 'figures/{2}_{0}_{1}_mean.pdf'.format(var,calendar.month_abbr[select_month],model)
    fig.savefig(outFile,transparent=True,bbox_inches='tight')
    print(outFile)
    # seasonal
    for select_season in select_seasons:
        tmp = dTS.sel(time=slice('{0:04d}-12'.format(avg_years[0]-1),'{0:04d}'.format(avg_years[1])))
        seas = tmp['time.season'] == select_season
        dTS_MAM = tmp[var].isel(time=seas).groupby('time.year').mean()
        # totals
        fig,ax,transf = ac.Projection('PlateCarree',kw_args={'central_longitude':180})
        dTS_tmp = dTS_MAM.mean('year')
        fig.set_figwidth(12)
        transf['cbar_kwargs'] = {'shrink':0.9,'label':labels[var]}
        pval = ac.StatTest(dTS_tmp,0,'T','member',parallel=True)
        if var in vmaxs.keys():
            transf['vmax'] = vmaxs[var]
        dTS_tmp.mean('member').plot(ax=ax,cmap=cmaps[var],**transf)
        del transf['cbar_kwargs']
        if var in vmaxs.keys():
            del transf['vmax']
        dTS_tmp.mean('member').where(pval<0.1).plot.contourf(ax=ax,colors='none',hatches=hatch,add_colorbar=False,**transf)
        ax.gridlines()
        ax.coastlines()
        ax.set_title('{0} anomalies, {1}, years {2}-{3}'.format(var,select_season,avg_years[0],avg_years[1]))
        outFile = 'figures/{2}_{0}_{1}_mean.pdf'.format(var,select_season,model)
        fig.savefig(outFile,transparent=True,bbox_inches='tight')
        print(outFile)
        # by year
        all_years = np.unique(dTS.time.dt.year)
        ncols = 5
        nrows = (len(all_years)-1)//ncols+1
        fig,axs,transf = ac.Projection('PlateCarree',ncols=ncols,nrows=nrows,kw_args={'central_longitude':180})
        transf['add_colorbar'] = False
        fig.set_figwidth(6*ncols)
        fig.set_figheight(3*nrows)
        for y,year in enumerate(all_years):
            ax = axs.flatten()[y]
            filtr = dTS.time.dt.season == select_season
            dTS_MAM = dTS[var].isel(time=filtr).groupby('time.year').mean()
            pval = ac.StatTest(dTS_MAM.sel(year=year),0,'T','member',parallel=True)
            cf = dTS_MAM.sel(year=year).mean('member').plot(ax=ax,vmin=mins[var],cmap=cmaps[var],**transf)
            dTS_MAM.sel(year=year).mean('member').where(pval<0.1).plot.contourf(ax=ax,colors='none',hatches=hatch,**transf)
            ax.gridlines()
            ax.coastlines()
            ax.set_title('{0} year {1}'.format(select_season,year))
        ac.AddColorbar(fig,axs,cf,shrink=0.5)
        outFile = 'figures/{2}_{0}_{1}.pdf'.format(var,select_season,model)
        fig.savefig(outFile,transparent=True,bbox_inches='tight')
        print(outFile)


def SaveFig(fig,name):
    fig.savefig(name,bbox_inches='tight',transparent=True)
    print(name)

def AddHatch(mam,mom,fig,axs,transf):
    #hatch=['//']
    #plt.rcParams['hatch.linewidth'] = 1.0
    #plt.rcParams['hatch.color'] = 'black'
    # plot hatches over mam where sign is the same as in mom
    months = mam.time.dt.month.values
    nmonth = len(months)
    if len(axs) != nmonth:
        print('NUMBER OF MONTHS != NUMBER OF AXES!')
        return
    for m,month in enumerate(months):
        ax = axs[m]
        mon = calendar.month_abbr[month]
        tmp = mam.isel(time=m).squeeze()
        tmp2= mom.isel(time=m).squeeze()
        # check where the sign is the same
        filtr = (tmp*tmp2) > 0
        ttle = ax.get_title()
        tmp.where(filtr).plot.contourf(ax=ax,colors='none',hatches=hatch,**transf)
        ax.set_title(ttle)

def PlotSeasonal(mam,season,fig_out=False):
    var = mam.name
    #if season == 'DJF':
    #    mam = mam.shift(time=3)
    #    mam = mam.sel(time=slice('0002',None))
    #    filtr = mam['time.season'] == 'MAM'
    #else:
    #    filtr = mam['time.season'] == season
    mam = mam.sel(time=slice('0001-12',None))
    filtr = mam['time.season'] == season
    mam = mam.isel(time=filtr).groupby('time.year').mean()
    nyears = len(mam.year)
    pval = ac.StatTest(mam,0,'T','member',parallel=True)
    mam = mam.mean('member')
    fig,axs,transf = ac.Projection('PlateCarree',ncols=nyears,kw_args={'central_longitude':180})
    fig.set_figwidth(6*nyears)
    for a,ax in enumerate(axs):
        cf=mam.isel(year=a).plot(ax=ax,cmap=cmaps[var],add_colorbar=False,vmin=mins[var],**transf)
        mam.where(pval<0.1).isel(year=a).plot.contourf(ax=ax,colors='none',add_colorbar=False,hatches=hatch,**transf)
        ax.gridlines()
        ax.coastlines()
        ax.set_title('Year {0}'.format(a+2))
    ac.AddColorbar(fig,axs,cf,shrink=0.3)
    if fig_out:
        return fig,axs,transf
    outFile = 'figures/{2}_{0}_{1}_yearly.pdf'.format(var,season,model)
    fig.savefig(outFile,bbox_inches='tight',transparent=True)
    print(outFile)


def PlotMonthly(mam,var,fig_out=False):
    months = mam.time.dt.month.values
    nmon = len(months)
    if nmon > 6:
        ncols = nmon//2
        nrows = (nmon-1)//ncols + 1
    else:
        ncols = nmon
        nrows = 1
    fig,axs,transf = ac.Projection('PlateCarree',ncols=ncols,nrows=nrows,kw_args={'central_longitude':180})
    fig.set_figwidth(6*ncols)
    transf['add_colorbar'] = False
    for m,month in enumerate(months):
        if nrows > 1:
            ax = axs.flatten()[m]
        else:
            ax = axs[m]
        # I don't use month  names because WACCM adds the perturbation on
        #  Jan 1 instead of Jan 15, so we want to use "Month 1" etc
        #mon = calendar.month_abbr[month]
        mon = 'month {0}'.format(month)
        tmp = mam.isel(time=m).squeeze()
        if 'member' in tmp.coords:
            pval = ac.StatTest(tmp,0,'T','member',parallel=True)
            tmp = tmp.mean('member')
        else:
            pval = xr.zeros_like(tmp)
        #transf['cbar_kwargs'] = {'shrink':0.5}
        cf = tmp.plot(ax=ax,cmap=cmaps[var],vmin=mins[var],center=0,**transf)
        #del transf['cbar_kwargs']
        if 'member' in mam.coords:
            tmp.where(pval<0.1).plot.contourf(ax=ax,colors='none',hatches=hatch,**transf)
        ax.gridlines()
        ax.coastlines()
        ax.set_title('{0}, {1} year 1'.format(var,mon))
    ac.AddColorbar(fig,axs,cf,shrink=0.5)
    if fig_out:
        return fig,axs,transf
    outFile = 'figures/{1}_{0}_maps.pdf'.format(var,model)
    SaveFig(fig,outFile)

waccm_slice = slice('0001-01','0001-03')
    
for var in dTS.data_vars:
    mam = dTS[var].sel(time=waccm_slice)
    PlotMonthly(mam,var)
#    for select_season in select_seasons:
#        PlotSeasonal(dTS[var],select_season)

# NOAA OLR
olr = xr.open_dataarray('/g/data/w40/mxj563/OLR/olr.mon.mean.nc')
olrm = xr.open_dataset('/g/data/w40/mxj563/OLR/olr.mon.ltm.1991-2020.nc').olr
olra = olr.sel(time='2022').groupby('time.month') - olrm.groupby('time.month').mean()

mam = olra.sel(time=slice('2022-01','2022-03'))
#PlotMonthly(mam,'OLR')
fig1,axs1,transf1 = PlotMonthly(mam,'OLR',fig_out=True)
mom = dTS[OLR].mean('member').sel(time=waccm_slice)
mom = mom.interp({'lon':mam.lon,'lat':mam.lat})
AddHatch(mam,mom,fig1,axs1,transf1)
SaveFig(fig1,'figures/OLR_maps.pdf')

# GPCP anomalies downloaded from WRIT
gpcp = xr.open_dataarray('gpcp_2022.nc')
#PlotMonthly(gpcp,'GPCP')
fig2,axs2,transf2 = PlotMonthly(gpcp,'GPCP',fig_out=True)
mom = dTS[PREC].mean('member').sel(time=waccm_slice)
mom = mom.interp({'lon':gpcp.lon,'lat':gpcp.lat})
AddHatch(gpcp,mom,fig2,axs2,transf2)
SaveFig(fig2,'figures/GPCP_maps.pdf')
