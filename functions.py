import xarray as xr
import pandas as pd
from aostools import climate as ac
import seaborn as sns
import numpy as np

def SaveFig(fig,outFile):
    fig.savefig(outFile,transparent=True,bbox_inches='tight')
    print(outFile)

def ReadFile(file):
    ds = xr.open_dataset(file,decode_times=False)
    attrs = ds.time.attrs
    ds = ds.assign_coords(time=(ds.time-4*365)/365)
    ds.time.attrs['units'] = 'years'
    #ds.time.attrs['calendar'] = attrs['calendar']
    ds = xr.decode_cf(ds)
    dt = ds.interp(time=1.0)
    ds = xr.concat([dt,ds],'time')
    return ds

def MakeDataFrame(ds,colname):
    memb = []
    years = []
    vals = []
    cls = []
    for col in ds[colname].values:
        dt = ds.sel({colname:col})
        for m in range(len(dt.member)):
            for y in range(len(dt.year)):
                cls.append(col)
                memb.append(dt.member.values[m])
                years.append(dt.year.values[y])
                vals.append(float(dt.isel(member=m,year=y).values))
    yr1 = years[0]
    years = [y-yr1+1 for y in years]
    df = pd.DataFrame(data={colname:cls,'member':memb,'year':years,ds.name:vals})
    return df

def CorrectTime(ds):
    calendar = ds.time.attrs['calendar']
    units = ds.time.attrs['units']
    if calendar == 'noleap': #waccm
        daysperyear = 365
        mid_month   = 0 #this is now fixed
    else:
        daysperyear = 360
        mid_month   = 0
    num_years = ds.time[0]//daysperyear
    ntime = ds.time - num_years*daysperyear + mid_month
    ntime.attrs = ds.time.attrs
    return xr.decode_cf(ds.assign_coords(time=ntime)),units,calendar


def GlobalMeanPlot(ta,name=None,fig_out=False,fig=None,ax=None,ttle=None,pval_parallel=True,plim=0.1,label='_none_',color=None):
    from matplotlib import pyplot as plt
    import cftime
    import nc_time_axis
    pval = ac.StatTest(ta,0,'T','member',parallel=pval_parallel)
#    pval = ac.StatTest(ta,None,'sign','member',parallel=pval_parallel)
    if fig is None and ax is None:
        fig,ax = plt.subplots()
        colr = 'r'
    else:
        colr = 'b'
    if color is not None:
        colr = color
    #ta.plot.line(ax=ax,x='time',color=colr,alpha=0.3,add_legend=False)
    ta_std = ta.std('member')
    ta_mn  = ta.mean('member')
    times = ta.time.values
    ax.plot(times,ta_mn,color=colr,ls='--',lw=1,label='_none_')
    if plim < 0.5:
        ta_mns = ta_mn.where(pval<plim)
    else:
        ta_mns = ta_mn.where(pval>plim)
    ax.plot(times,ta_mns,color=colr,lw=2,label=label)
    #ta_mn.plot.line(ax=ax,x='time',color=colr,ls='--',lw=1)
    #ta_mn.where(pval<0.1).plot.line(ax=ax,x='time',color=colr,lw=2)
    lower = ta_mn - ta_std
    upper = ta_mn + ta_std
    ax.fill_between(times,lower,upper,color=colr,alpha=0.3,label='_none_')
    first_month = ta.time.dt.month.values[0]
    first_tick = (6-first_month+1)%6
    tick_vals = times[first_tick::6]
    ax.set_xticks(tick_vals)
    #xlims = ax.get_xlim()
    #if first_month == 1:
    #    xlims[0] = times[0]
    #else:
    #    xlims[0] = cftime.datetime(times[0].year,1,times[0].day)
    #ax.set_xlim(cftime.datetime(times[0].year,1,1),xlims[1])
    labls = []
    for tick in tick_vals:
        if tick.month == 7:
            labls.append('')
        else:
            labls.append('{}'.format(tick.year-1))
    ax.set_xticklabels(labls)
    ax.set_xlabel('time [years since eruption]')
    if ttle is None:
        ax.set_title('Global mean surface temperature')
    else:
        ax.set_title(ttle)
    ax.axhline(0,color='k')
    ax.grid()
    sns.despine()
    if fig_out:
        return fig,ax
    if name is not None:
        outFile = 'figures/{0}.pdf'.format(name)
        fig.savefig(outFile,bbox_inches='tight',transparent=True)
        print(outFile)

def ConvertTime2Years(ds):
    fyr = []
    for t in ds.time:
        fyr.append(t.dt.year-1 + t.dt.dayofyear/365)
    fyr = xr.DataArray(fyr,[('time',fyr)],name='time')
    fyr.attrs['units'] = 'years'
    return ds.assign_coords(time=fyr)


def ReadMLS(pure_anom=False):
    mls  = xr.open_dataset('MLS_data.nc')
    mls = mls.assign_coords({'TS_Time':mls.ts_Time,'ZM_Lat':mls.zm_lat})
    mls = mls.rename({'TS_Time':'time','ZM_Lat':'lat'})
    # align real eruption (max around Jan 25) with WACCM sims (Jan 1)
    mls.time.attrs['units'] = 'days since 0000-10-05'
    mls.time.attrs['calendar'] = 'noleap'
    mls = xr.decode_cf(mls)
    if pure_anom:
        base = mls.sel(time='0000-12-19').mean('time')
        return mls - base
    else:
        return mls

def ReadMLSMap():
    mls = xr.open_dataset('MLS_data.nc')
    time_name = None
    lon_name = None
    lat_name = None
    renames = {}
    for coord in mls.coords:
        if 'time' in coord.lower():
            time_name = coord
            renames[coord] = 'time'
        elif 'lon' in coord.lower():
            mls = mls.isel({coord:slice(0,-1)})
            lon_name = coord
            renames[coord] = 'lon'
        elif 'lat' in coord.lower():
            lat_name = coord
            renames[coord] = 'lat'
    #mls = mls.assign_coords({'Map_Time':mls.map_time,'Map_Long':mls.map_lon,'Map_Lat':mls.map_lat})
    mls = mls.rename(renames)
    mls.time.attrs['units'] = 'days since 0001-01-01'
    mls.time.attrs['calendar'] = 'noleap'
    mls = xr.decode_cf(mls)
    return mls.map_data


variables = {'waccm': {'T':'T','U':'U','SLP':'PSL','O3':'O3','TS':'TS','TREFHT':'TREFHT','Q':'Q','P':'PREC','TCWV':'TMQ','TCO':'TCO','OLR':'FLNT','OMEGA':'OMEGA','TH':'TH','VTH3d':'VTH3d'},
             'mima' : {'T':'temp','U':'ucomp','SLP':'slp','TS':'t_surf','Q':'sphum','P':'precip','OLR':'olr','TCWV':'tcwv'}}
for mod in ['aqua','aqua_sponge','aqua_sponge_10yr','bench_SH']:
    variables[mod] = variables['mima']


def Mass(ds):
    '''mass above a given pressure level, in kg/m2.
    '''
    from aostools.constants import g,a0,coslat
    mass_y = ds
    #mass_y = coslat(ds['lat'])*ds
    #mass_y = a0*ds/np.deg2rad(coslat(ds['lat']).integrate('lat'))
    #mass_p = mass_y.mean('lon').integrate('lev')*100
    mass_p = mass_y.integrate('lev')*100
    mass = mass_p/g
    if ds.lev[0] > ds.lev[-1]:
        mass = -mass
    return mass

def ExpFit(x,y):
    '''y= A*exp(B*x), x is assumed in units of years.
    '''
    fit = np.polyfit(x,np.log(y),1)
    lifetime = -1/fit[0]
    halftime = np.log(0.5)/fit[0]
    percent_per_month = (1-np.exp(fit[0]*1/12))*100
    percent_per_year  = (1-np.exp(fit[0]*1   ))*100
    return lifetime,halftime,percent_per_month,percent_per_year,fit[1]


def AddBox(ax,region,colr):
    from shapely import geometry
    from cartopy import crs as ccrs
    if region['lon'].start > 180:
        lonstart = region['lon'].start - 360
    else:
        lonstart = region['lon'].start
    if region['lon'].stop > 180:
        lonstop = region['lon'].stop - 360
    else:
        lonstop = region['lon'].stop
    geom = geometry.box(minx=lonstart,maxx=lonstop,
                                    miny=region['lat'].start,maxy=region['lat'].stop)
    #geoms.append(geom)
    ax.add_geometries([geom],edgecolor=colr,facecolor='none',crs=ccrs.PlateCarree(),zorder=5)
