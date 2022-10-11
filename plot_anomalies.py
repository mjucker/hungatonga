import xarray as xr
from aostools import climate as ac
from matplotlib import pyplot as plt
import seaborn as sns
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',help='Choose model/run to plot.')
parser.add_argument('-t',dest='topo',action='store_true',help='Whether or not to add continents in figures.')
parser.add_argument('-v',dest='vars',nargs='+',default=None,help='Limit to these variables.')
args = parser.parse_args()

model = args.model

pval_parallel = True

variables = {'waccm': {'T':'T','U':'U','SLP':'PSL','O3':'O3','TS':'TS','Q':'Q','P':'PREC','TCWV':'TMQ'},
             'mima' : {'T':'temp','U':'ucomp','SLP':'slp','TS':'t_surf','Q':'sphum','P':'precip','OLR':'olr','TCWV':'tcwv'}}

limits = {'Q': {'pres':slice(100,1)}}

nlevs = 16

# assign MiMA experiments to MiMA
for cse in ['aqua_sponge_5yr','bench_SH']:
    variables[cse] = variables['mima']

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

delta = DiscardVars(xr.open_dataset(pert_file)) - DiscardVars(xr.open_dataset(ctrl_file))
delta = ac.StandardGrid(delta,rename=True,pdir='decrease')
#if delta.pres[0] < delta.pres[-1]:
#    delta = delta.isel(pres=(None,None,-1))

def AddDates(ds):
    #from datetime import timedelta
    #WACCM data starts on February
    ## add NaNs to January of the first year
    if ds.time.dt.month.values[0] == 2:
        yr = ds.time.dt.year.values[0]
        #nans = ds.isel(time=0)/0
        #nans.assign_coords(time=nans.time-timedelta(days=31))
        #print(nans)
        #ds = xr.concat([nans,ds],'time')
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
            ax.set_title('{0}, year {1:d}'.format(calendar.month_abbr[mon],year))

def FixLabels(fig):
    from matplotlib.ticker import ScalarFormatter
    axs = fig.axes
    nrows = len(axs)
    ncols = len(axs[0])
    for i in range(nrows):
        for j in range(ncols):
            if j>0:
                axs[i][j].set_ylabel('')
                #axs[i][j].set_yticklabels('')
            if i<nrows-1:
                axs[i][j].set_xlabel('')
                #axs[i][j].set_xticklabels('')
    axs[0][0].get_yaxis().set_major_formatter(ScalarFormatter())
        
    
hatch = ['..']

for var in plot_vars:
     if var in plot_limits.keys():
         tmp = delta[var].sel(plot_limits[var])
     else:
         tmp = delta[var]
     tmp = AddDates(tmp)
     if 'pres' in tmp.coords: 
         tmp = tmp.sel(pres=slice(None,1)).mean('lon')
         pval = ac.StatTest(tmp,0,'T','member',parallel=pval_parallel) 
         fig = tmp.mean('member').plot.contourf(levels=nlevs,x='lat',col='time',col_wrap=12,yscale='log',yincrease=False,robust=True)
         for a,t in enumerate(tmp.time):
             ax = fig.axes.flatten()[a]
             tmp.mean('member').where(pval>0.1).isel(time=a).plot.contourf(ax=ax,x='lat',yscale='log',yincrease=False,hatches=hatch,colors='none',add_colorbar=False)
     else:
         pval = ac.StatTest(tmp,0,'T','member',parallel=pval_parallel)
         fig = tmp.mean('member').plot.contourf(levels=nlevs,x='lon',col='time',col_wrap=12,robust=True) 
         for a,t in enumerate(tmp.time):
             ax = fig.axes.flatten()[a]
             tmp.mean('member').where(pval>0.1).isel(time=a).plot.contourf(ax=ax,x='lon',hatches=hatch,colors='none',add_colorbar=False)
             if args.topo:
                 topo.plot.contour(levels=[0],colors='k',x='lon',add_colorbar=False,ax=ax)
             else:
                 ax.grid()
     ChangeTitles(fig)
     FixLabels(fig)
     ivar = inv_vars[var]
     outFile = 'figures/{0}_{1}_monthly.pdf'.format(model,ivar) 
     plt.savefig(outFile,bbox_inches='tight',transparent=True) 
     print(outFile) 


def GlobalMeanPlot(ta,name=None,fig_out=False,fig=None,ax=None):
    pval = ac.StatTest(ta,0,'T','member',parallel=pval_parallel)
    if fig is None and ax is None:
        fig,ax = plt.subplots()
        colr = 'r'
    else:
        colr = 'b'
    #ta.plot.line(ax=ax,x='time',color=colr,alpha=0.3,add_legend=False)
    ta_std = ta.std('member')
    ta_mn  = ta.mean('member')
    ta_mn.plot.line(ax=ax,x='time',color=colr,ls='--',lw=1)
    ta_mn.where(pval<0.1).plot.line(ax=ax,x='time',color=colr,lw=2)
    #(ta_mn-ta_std).plot.line(ax=ax,x='time',color=colr,ls='--')
    #(ta_mn+ta_std).plot.line(ax=ax,x='time',color=colr,ls='--')
    lower = ta_mn - ta_std
    upper = ta_mn + ta_std
    xdata = ax.get_lines()[0].get_xdata() 
    plt.fill_between(xdata,lower,upper,color=colr,alpha=0.3)
    ax.set_title('Global mean surface temperature')
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
    ta = ac.GlobalAvgXr(delta[tsurf].mean('lon'),[-90,90])
    GlobalMeanPlot(ta,'{0}_TS_global'.format(model))
    ta = ac.GlobalAvgXr(delta[tsurf].mean('lon'),[-90,0])
    tb = ac.GlobalAvgXr(delta[tsurf].mean('lon'),[0,90])
    fig,ax = GlobalMeanPlot(ta,fig_out=True)
    GlobalMeanPlot(tb,name=None,fig=fig,ax=ax)
    ax.grid()
    ax.text(0.05,0.90,'-- SH',color='r',transform=ax.transAxes)
    ax.text(0.05,0.85,'-- NH',color='b',transform=ax.transAxes)
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
    outFile = 'figures/{0}_AnnularModes.pdf'.format(model)
    fig.savefig(outFile,bbox_inches='tight',transparent=True)
    print(outFile)
    
