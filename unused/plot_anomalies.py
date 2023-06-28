import xarray as xr
from aostools import climate as ac
from matplotlib import pyplot as plt
from hungatonga import functions as fc
import seaborn as sns
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',help='Choose model/run to plot.')
parser.add_argument('-t',dest='topo',action='store_true',help='Whether or not to add continents in figures.')
parser.add_argument('-v',dest='vars',nargs='+',default=None,help='Limit to these variables.')
parser.add_argument('-Y',dest='years',nargs=2,type=int,default=[None,None],help='Limit to these years.')
parser.add_argument('-M',dest='months',nargs=2,type=int,default=[None,None],help='Limit to these months.')
parser.add_argument('-r',dest='robust',action='store_true',help='Scale colorbar robustly or not.')
parser.add_argument('-n',dest='num',type=int,default=16,help='Number of color levels to use.')
parser.add_argument('-l',dest='lat_lim',type=float,nargs=2,default=[None,None],help='Limit to these latitudes.')
args = parser.parse_args()

sns.set_context('paper')#,font_scale=1.5)

model = args.model
years = [None,None]
for y,yr in enumerate(args.years):
    if yr is not None:
        years[y] = '{0:04d}'.format(yr)
if years[0] is None and years[1] is None:
    for m,mnth in enumerate(args.months):
        if mnth is not None:
            years[m] = '0001-{0:02d}'.format(mnth)


pval_parallel = True

variables = fc.variables #{'waccm': {'T':'T','U':'U','SLP':'PSL','O3':'O3','TS':'TS','Q':'Q','P':'PREC','TCWV':'TMQ','TCO':'TCO','OLR':'FLNT'},
             #'mima' : {'T':'temp','U':'ucomp','SLP':'slp','TS':'t_surf','Q':'sphum','P':'precip','OLR':'olr','TCWV':'tcwv'}}            

limits = {'Q': {'pres':slice(100,1)}}
vmins  = {'Q': 0}

scale = {'Q' : 1e6, 'P': 1e3*86400}


cmaps = {
    'TS': 'RdBu_r',
    'P' : 'PuOr',
    'OLR':'PiYG',
    'SLP':'BrBG_r',
    'Q'  :'Blues',
    'U'  : 'RdYlBu_r',
    }
    


nlevs = args.num
robust = args.robust

# assign MiMA experiments to MiMA
for cse in ['aqua_sponge_5yr','bench_SH','aqua_sponge_10yr','bench_SH_oz']:
    variables[cse] = variables['mima']
# assign WACCM experiments to WACCM
for cse in ['waccm_PI','waccm_F2000','waccm_F2000_10yrs']:
    variables[cse] = variables['waccm']

base = '{0}_xxx_ens.nc'.format(model)
ctrl_file = base.replace('xxx','ctrl')
pert_file = base.replace('xxx','pert')

if args.vars is None:
    plot_vars = variables[args.model].values()
    plot_limits = {}
    for var in variables[args.model].keys():
        if var in limits.keys():
            plot_limits[variables[args.model][var]] = limits[var]
else:
    plot_vars = [variables[args.model][v] for v in args.vars]
    plot_limits = {}
    for var in args.vars:
        if var in limits.keys():
            plot_limits[variables[args.model][var]] = limits[var]
print('PLOTTING THE FOLLOWING VARIABLE(S): ',plot_vars)
inv_vars = {v:k for k,v in variables[args.model].items()}

def DiscardVars(ds):
    disc_vars = []
    for var in ds.data_vars:
        if var not in plot_vars:
            disc_vars.append(var)
    for var in disc_vars:
        del ds[var]
    return ds

#def CorrectTime(ds):
#    calendar = ds.time.attrs['calendar']
#    units = ds.time.attrs['units']
#    if calendar == 'noleap': #waccm
#        daysperyear = 365
#        #mid_month   = -15 # this is now fixed in analyze_ht_monthly_waccm.py
#        mid_month   = 0
#    else:
#        daysperyear = 360
#        mid_month   = 0
#    num_years = ds.time[0]//daysperyear
#    ntime = ds.time - num_years*daysperyear + mid_month
#    ntime.attrs = ds.time.attrs
#    return xr.decode_cf(ds.assign_coords(time=ntime)),units,calendar

delta = DiscardVars(xr.open_dataset(pert_file,decode_times=False)) - DiscardVars(xr.open_dataset(ctrl_file,decode_times=False))
delta = ac.StandardGrid(delta,rename=True,pdir='decrease')
delta = delta.sel(lat=slice(*args.lat_lim))

# total column ozone
if 'TCO' not in delta.data_vars and 'TCO' in plot_vars:
    if 'O3' in variables[args.model].keys():
        ozone = variables[args.model]['O3']
    else:
        ozone = None
    if 'T' in  variables[args.model].keys():
        temp  = variables[args.model]['T']
    else:
        temp = None
    if ozone in plot_vars and temp in plot_vars:
        print('COMPUTING TOTAL COLUMN OZONE')
        o3_ctrl = xr.open_dataset(ctrl_file,decode_times=False)[ozone]
        T_ctrl  = xr.open_dataset(ctrl_file,decode_times=False)[temp]
        o3_pert = xr.open_dataset(pert_file,decode_times=False)[ozone]
        T_pert  = xr.open_dataset(pert_file,decode_times=False)[temp]
        tco_ctrl= ac.TotalColumnOzone(o3_ctrl,T_ctrl)
        tco_pert= ac.TotalColumnOzone(o3_pert,T_pert)
        delta['TCO'] = tco_pert-tco_ctrl
delta,units,calendar = fc.CorrectTime(delta)

#if delta.pres[0] < delta.pres[-1]:
#    delta = delta.isel(pres=(None,None,-1))

def AddDates(ds):
    #from datetime import timedelta
    #WACCM data starts on February
    ## add NaNs to January of the first year
    if ds.time.dt.month.values[0] == 2:
        yr = ds.time.dt.year.values[0]
        ds = xr.concat([ds.interp(time='{0:04d}-01-01'.format(yr)),ds],'time')
    return ds


if args.topo:
    topo = xr.open_dataarray('zsurf.nc')
    topo = topo.interp({'lon':delta.lon,'lat':delta.lat})

def ChangeTitles(fig):
    import calendar
    axs = fig.axes
    first_year = None
    for ax in axs.flatten():
        ttle = ax.get_title()
        ttle_split = ttle.split('-')
        if len(ttle_split)>1:
            mon = int(ttle_split[1])
            year= int(ttle.split('=')[1].split('-')[0])
            if first_year is None:
                first_year = year
            year = year-first_year+1
            ax.set_title('month {0}, year {1:d}'.format(mon,year))

def FixLabels(fig):
    from matplotlib.ticker import ScalarFormatter
    axs = fig.axes
    nrows = len(axs)
    ncols = len(axs[0])
    for i in range(nrows):
        for j in range(ncols):
            if 'lat' in axs[i][j].get_xlabel():
                axs[i][j].axvline(0,color='k',lw=1)
            if j>0:
                axs[i][j].set_ylabel('')
                #axs[i][j].set_yticklabels('')
            if i<nrows-1:
                axs[i][j].set_xlabel('')
                #axs[i][j].set_xticklabels('')
            else:
                labl = axs[i][j].get_xlabel()
                if 'lat' in labl and '[' in labl:
                   axs[i][j].set_xlabel('lat')
    axs[0][0].get_yaxis().set_major_formatter(ScalarFormatter())
        
    
hatch = ['..']
plim = 0.1
for var in plot_vars:
     if inv_vars[var] in plot_limits.keys():
         tmp = delta[var].sel(plot_limits[var])
     else:
         tmp = delta[var]
     if inv_vars[var] in scale.keys():
         tmp = scale[inv_vars[var]]*tmp
     if inv_vars[var] in cmaps.keys():
         cmap = cmaps[inv_vars[var]]
     else:
         cmap = None
     if inv_vars[var] in vmins.keys():
         vmin = vmins[inv_vars[var]]
     else:
         vmin = None
     #tmp = AddDates(tmp)
     tmp = tmp.sel(time=slice(*years))
     if 'pres' in tmp.coords: 
         tmp = tmp.sel(pres=slice(None,1)).mean('lon')
         pval = ac.StatTest(tmp,0,'T','member',parallel=pval_parallel)
#         pval = ac.StatTest(tmp,None,'sign','member',parallel=pval_parallel)
         fig = tmp.mean('member').plot.contourf(levels=nlevs,x='lat',col='time',col_wrap=min(12,len(tmp.time)),yscale='log',yincrease=False,cmap=cmap,vmin=vmin,robust=robust)
         for a,t in enumerate(tmp.time):
             ax = fig.axes.flatten()[a]
             if plim < 0.5:
                 tmps = tmp.mean('member').where(pval>plim)
             else:
                 tmps = tmp.mean('member').where(pval<plim)
             tmps.isel(time=a).plot.contourf(ax=ax,x='lat',yscale='log',yincrease=False,hatches=hatch,colors='none',add_colorbar=False)
     else:
         pval = ac.StatTest(tmp,0,'T','member',parallel=pval_parallel)
#         pval = ac.StatTest(tmp,None,'sign','member',parallel=pval_parallel)
         fig = tmp.mean('member').plot.contourf(levels=nlevs,x='lon',col='time',col_wrap=min(12,len(tmp.time)),cmap=cmap,vmin=vmin,robust=robust) 
         for a,t in enumerate(tmp.time):
             ax = fig.axes.flatten()[a]
             if plim < 0.5:
                 tmps = tmp.mean('member').where(pval>plim)
             else:
                 tmps = tmp.mean('member').where(pval<plim)
             tmps.isel(time=a).plot.contourf(ax=ax,x='lon',hatches=hatch,colors='none',add_colorbar=False)
             if args.topo:
                 topo.plot.contour(levels=[0],colors='k',x='lon',add_colorbar=False,ax=ax)
             else:
                 ax.grid()
     #plt.tight_layout()
     ChangeTitles(fig)
     FixLabels(fig)
     ivar = inv_vars[var]
     outFile = 'figures/{0}_{1}_monthly.pdf'.format(model,ivar) 
     plt.savefig(outFile,transparent=True,bbox_inches='tight') 
     print(outFile) 


def GlobalMeanPlot(ta,name=None,fig_out=False,fig=None,ax=None,ttle=None):
    import cftime
    import nc_time_axis
    pval = ac.StatTest(ta,0,'T','member',parallel=pval_parallel)
#    pval = ac.StatTest(ta,None,'sign','member',parallel=pval_parallel)
    if fig is None and ax is None:
        fig,ax = plt.subplots()
        colr = 'r'
    else:
        colr = 'b'
    #ta.plot.line(ax=ax,x='time',color=colr,alpha=0.3,add_legend=False)
    ta_std = ta.std('member')
    ta_mn  = ta.mean('member')
    times = ta.time.values
    ax.plot(times,ta_mn,color=colr,ls='--',lw=1)
    if plim < 0.5:
        ta_mns = ta_mn.where(pval<plim)
    else:
        ta_mns = ta_mn.where(pval>plim)
    ax.plot(times,ta_mns,color=colr,lw=2)
    #ta_mn.plot.line(ax=ax,x='time',color=colr,ls='--',lw=1)
    #ta_mn.where(pval<0.1).plot.line(ax=ax,x='time',color=colr,lw=2)
    lower = ta_mn - ta_std
    upper = ta_mn + ta_std
    plt.fill_between(times,lower,upper,color=colr,alpha=0.3)
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
    ax.set_xlabel('time [years]')
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


tsurf = variables[model]['TS']
if tsurf in delta.data_vars:
    #waccm has fixed SSTs, so concentrate on land only
    if 'waccm' in model.lower():
        mask = xr.open_dataarray('landmask_waccm.nc')
        delta[tsurf] = delta[tsurf]*mask
    ta = ac.GlobalAvgXr(delta[tsurf].mean('lon'),[-90,90])
    GlobalMeanPlot(ta,'{0}_TS_global'.format(model))
    ta = ac.GlobalAvgXr(delta[tsurf].mean('lon'),[-90,-30])
    tb = ac.GlobalAvgXr(delta[tsurf].mean('lon'),[30,90])
    fig,ax = GlobalMeanPlot(ta,fig_out=True)
    GlobalMeanPlot(tb,name=None,fig=fig,ax=ax)
    ax.grid()
    ax.text(0.05,0.90,'-- SH 30-90',color='r',transform=ax.transAxes)
    ax.text(0.05,0.85,'-- NH 30-90',color='b',transform=ax.transAxes)
    outFile = 'figures/{0}_TS_hemi.pdf'.format(model)
    fig.savefig(outFile,bbox_inches='tight',transparent=True)
    print(outFile)
    #
slp = variables[model]['SLP']
if slp in delta.data_vars:
    sam = delta[slp].mean('lon').sel(lat=-40,method='nearest') \
        - delta[slp].mean('lon').sel(lat=-60,method='nearest')
    fig,ax = GlobalMeanPlot(sam,fig_out=True)
    nam = delta[slp].mean('lon').sel(lat=40,method='nearest') \
        - delta[slp].mean('lon').sel(lat=60,method='nearest')
    GlobalMeanPlot(nam,name=None,fig=fig,ax=ax)
    ax.grid()
    ax.text(0.05,0.90,'-- SAM',color='r',transform=ax.transAxes)
    ax.text(0.05,0.85,'-- NAM',color='b',transform=ax.transAxes)
    ax.set_title('Annular Modes')
    outFile = 'figures/{0}_AnnularModes.pdf'.format(model)
    fig.savefig(outFile,bbox_inches='tight',transparent=True)
    print(outFile)

if 'TCO' in delta.data_vars and 'TCO' in plot_vars:
    nh = [80,90]
    sh = [-70,-50]
    tco = delta['TCO'].mean('lon')
    fig,ax = GlobalMeanPlot(ac.GlobalAvgXr(tco,sh),fig_out=True)
    #fig,ax = GlobalMeanPlot(tco.sel(lat=slice(-90,-50)).sum('lat'),fig_out=True)
    GlobalMeanPlot(ac.GlobalAvgXr(tco,nh),name=None,fig=fig,ax=ax)
    #GlobalMeanPlot(tco.sel(lat=slice(-90,-50)).sum('lat'),name=None,fig=fig,ax=ax)
    ax.grid()
    ax.text(0.05,0.90,'-- SH {1}-{0}'.format(abs(sh[0]),abs(sh[1])),color='r',transform=ax.transAxes)
    ax.text(0.05,0.85,'-- NH {0}-{1}'.format(*nh),color='b',transform=ax.transAxes)
    ax.set_title('Total Column Ozone [DU]')
    outFile = 'figures/{0}_TCO_hemi.pdf'.format(model)
    fig.savefig(outFile,bbox_inches='tight',transparent=True)
    print(outFile)

#def GlobalMass(da):
#    from aostools import constants as at
#    mass_y = ac.GlobalAvgXr(da,[-90,90])
#    mass_x = np.deg2rad(mass_y.integrate('lon'))
#    mass_p = mass_x.integrate('pres')*100
#    mass = at.a0**2/at.g*mass_p
#    if da['pres'][0]>da['pres'][-1]:
#        mass = -mass
#    return mass

#sphum = variables[model]['Q']
#if sphum in delta.data_vars:
#    sphum = delta[sphum]
#    if 'Q' in limits.keys():
#        sphum = sphum.sel(limits['Q'])
#    tq = ac.GlobalMass(sphum)
#    GlobalMeanPlot(tq*1e-9,'{0}_tQ'.format(model),ttle='Total stratospheric water mass [Tg]')

