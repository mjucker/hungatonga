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
    #yr1 = years[0]
    #years = [y-yr1+1 for y in years]
    df = pd.DataFrame(data={colname:cls,'member':memb,'year':years,ds.name:vals})
    return df

def CorrectTime(ds,decode=True):
    calendar = ds.time.attrs['calendar']
    units = ds.time.attrs['units']
    if calendar.lower() == 'noleap': #waccm
        daysperyear = 365
        mid_month   = 0 #this is now fixed
    else:
        daysperyear = 360
        mid_month   = 0
    if ds.time[0] <= 0:
        startyear = int(units.split('since ')[-1].split('-')[0])
        startmonth= int(units.split('-')[1].split('-')[0])
        if startmonth == 12: #start with DJF
            ds = ds.assign_coords({'time':ds.time+31})
            ds.time.attrs['calendar'] = calendar
        #    correctyear = 0
        #else:
        #    correctyear = 1
        correctyear = 0
        units = units.replace('{:04d}'.format(startyear),'{:04d}'.format(correctyear))
        ds.time.attrs['units'] = units
        if decode:
            ds = xr.decode_cf(ds)
        return ds,units,calendar
    else:
        num_years = ds.time[0]//daysperyear
        ntime = ds.time - num_years*daysperyear + mid_month
        ntime.attrs = ds.time.attrs
        ds = ds.assign_coords(time=ntime)
    if decode:
        ds = xr.decode_cf(ds)
    return ds,units,calendar


def GlobalMeanPlot(ta,name=None,fig_out=False,fig=None,ax=None,ttle=None,pval_parallel=True,plim=0.1,label='_none_',color=None,ls=None,fill=True):
    from matplotlib import pyplot as plt
    import cftime
    import nc_time_axis
    ta = ta.convert_calendar('noleap',align_on='date')
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
    ax.plot(times,ta_mns,color=colr,lw=2,label=label,ls=ls)
    #ta_mn.plot.line(ax=ax,x='time',color=colr,ls='--',lw=1)
    #ta_mn.where(pval<0.1).plot.line(ax=ax,x='time',color=colr,lw=2)
    lower = ta_mn - ta_std
    upper = ta_mn + ta_std
    if fill:
        ax.fill_between(times,lower,upper,color=colr,alpha=0.3,label='_none_')
    else:
        ax.plot(times, lower, color=colr,lw=0.5,alpha=0.7,label='_none_')
        ax.plot(times, upper, color=colr,lw=0.5,alpha=0.7,label='_none_')
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
            labls.append('{}'.format(tick.year))
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
    #mls = mls.assign_coords({'TS_Time':mls.ts_Time,'ZM_Lat':mls.zm_lat})
    if 'ts_Time' in mls.coords:
        mls = mls.rename({'ts_Time':'TS_Time'})
    if 'zm_lat' in mls.data_vars:
        mls = mls.rename({'zm_lat':'ZM_Lat'})
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


variables = {'waccm': {'T':'T','U':'U','V':'V','Z':'Z3','SLP':'PSL','O3':'O3','TS':'TREFHT','TREFHT':'TREFHT','Q':'Q','P':'PREC','TCWV':'TMQ','TCO':'TCO','OLR':'FLNT','DLS':'FLDS','DSS':'FSDS','OMEGA':'OMEGA','TH':'TH','VTH3d':'VTH3d','CLDTOT':'CLDTOT','CLDLOW':'CLDLOW','CLDHGH':'CLDHGH','CLDMED':'CLDMED','LWCF':'LWCF','SWCF':'SWCF','ICEFRAC':'ICEFRAC'},
             'mima' : {'T':'temp','U':'ucomp','V':'vcomp','SLP':'slp','TS':'t_surf','Q':'sphum','P':'precip','OLR':'olr','TCWV':'tcwv','DLS':'flux_lw','DSS':'flux_sw'}}
for mod in ['aqua','aqua_sponge','hthh','hthh_fix','bench_SH']:
    variables[mod] = variables['mima']


def ModName(model):
    if '_' in model:
        splitname = model.split('_')
        if len(splitname) > 2:
            qmodel = '_'.join(splitname[:2])
        else:
            qmodel = '_'.join(splitname[:-1])
    else:
        qmodel = model
    return qmodel
    
def Mass(ds):
    '''mass above a given pressure level, in kg/m2.
    '''
    from aostools.constants import g,a0,coslat
    from aostools import climate as ac
    names = ac.FindCoordNames(ds)
    lev = names['pres']
    mass_y = ds
    #mass_y = coslat(ds['lat'])*ds
    #mass_y = a0*ds/np.deg2rad(coslat(ds['lat']).integrate('lat'))
    #mass_p = mass_y.mean('lon').integrate('lev')*100
    mass_p = mass_y.integrate(lev)*100
    mass = mass_p/g
    if ds[lev][0] > ds[lev][-1]:
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

def CheckQBO(ens_file,model):
    from aostools import climate as ac
    U = variables[model]['U']
    if isinstance(ens_file,str):
        ds = xr.open_dataset(ens_file)[U]
    else:
        ds = ens_file[U]
    lev= ac.FindCoordNames(ds)['pres']
    U50 = ds.sel({lev:50,'lat':slice(-5,5)}).isel(time=0).mean(['lon','lat'])
    del U50['time']
    del U50[lev]
    qbo_pos = U50 > 0
    qbo_neg = U50 < 0
    return qbo_pos,qbo_neg

def RenameQBOFile(filename,qbo,ftype='.pdf'):
    qnme = {'+':'p','-':'m'}
    return filename.replace(ftype,'_QBO{0}{1}'.format(qnme[qbo],ftype))


areas = {}
areas['P']  = {'DJF':{
                #'ITCZ':{'lon':slice(180.1,210), 'lat': slice(5,10)},
                'Pacific':{'lon':slice(180.1,220),'lat':slice(5,25)},
                #'MC'  :{'lon':slice(120,140),   'lat': slice(0,20)},
                'MC'  :{'lon':slice(110,150),   'lat': slice(0,20)},
                #'Europe':{'lon':slice(0.1,30),   'lat': slice(45,53)},
                },
               'JJA':{
                #'ITCZ':{'lon':slice(180,210), 'lat': slice(5,10)},
                #'MC'  :{'lon':slice(120,140),   'lat': slice(0,20)},
                'ION':  {'lon':slice(65,95),     'lat': slice(10,20)},
                'IOS':  {'lon':slice(40,100),     'lat': slice(-10,-2)},
                },
               'MAM':{
                #'ITCZS':{'lon':slice(180.1,210), 'lat': slice(-15,-10)},
                #'SPCZN':{'lon':slice(145,165), 'lat': slice(-15,-10)},
                },
               'SON':{
                #'ITCZ':{'lon':slice(180.1,200), 'lat': slice(10,15)},
                #'MC'  :{'lon':slice(120,140),   'lat': slice(0,20)},
                },
               }
areas['TS'] = {'DJF':{
                'Scandinavia':  {'lon':slice(10,40)  ,'lat':slice(58,70)},
                'Eurasia': {'lon':slice(40,80)  ,'lat':slice(35,50)},
                'NAmerica':{'lon':slice(235,265),'lat':slice(45,65)},
                #'Australia':{'lon':slice(120,145),'lat':slice(-28,-18)},
                #'NAfrica': {'lon':slice(-10,30)  ,'lat':slice(20,35)},
               },
         'JJA':{
                #'Scandinavia':  {'lon':slice(20,60)  ,'lat':slice(55,70)},
                #'NAmerica': {'lon':slice(260,290),'lat':slice(40,60)},
                'Australia':{'lon':slice(120,145),'lat':slice(-28,-18)},
                'Arctic':{'lon':slice(40,270),'lat':slice(75,90)},
                'Amundsen':{'lon':slice(185,280),'lat':slice(-80,-65)},
               },
         'MAM':{
                'Siberia': {'lon':slice(50,90),'lat':slice(45,65)},
                'Arctic':{'lon':slice(40,270),'lat':slice(75,90)},
               },
         'SON': {
                'NAsia'  : {'lon':slice(70,120),'lat':slice(45,55)},
                'Arctic':{'lon':slice(40,270),'lat':slice(75,90)},
               },
         #'Arctic' : {'lon':slice(0,360),  'lat':slice(80,90)},
         #'EAsia':   {'lon':slice(100,125),'lat':slice(40,55)},
         #'Australia':{'lon':slice(115,155),'lat':slice(-38,-20)}
         }
for var in ['CLDTOT','DLS','OLR','LWCF','SWCF']:
    areas[var] = areas['TS']
