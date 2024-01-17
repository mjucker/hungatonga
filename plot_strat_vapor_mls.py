#!env python

import xarray as xr
from aostools import climate as ac
from aostools import constants as at
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from hungatonga import functions as fc
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-L',dest='level',default=None,help='Plot integrated water vapor [None], water vapor at a given level [level in hPa], or zonal mean of a given month [YYYY-MM], where YYYY is in MLS time, starting with 0001.')
args = parser.parse_args()


sns.set_context('notebook')
sns.set_style('whitegrid')
colrs = sns.color_palette('bright')

def MassPerDeg(ds):
    from aostools.constants import a0,g,coslat
    mass_y = coslat(ds.lat)*ds
    mass_x = np.deg2rad(mass_y.integrate('lon'))
    mass_p = mass_x.integrate(lev)*100
    mass = a0**2/g*mass_p
    if ds[lev][0] > ds[lev][-1]:
        mass = -mass
    return mass

def PlotMass(mass,tco=None,appendix='',levels=20,colr=None,fig_out=False,kind='auto',fig=None,ax=None,do_pval=True,ttle=None):
    #fig,ax = plt.subplots()
    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = fig.add_axes([0.1,0.1,0.74,0.8])
    fig.set_figwidth(8)
    cmass = None
    qmass = mass
    units = qmass.attrs['units']
    pmass = qmass
    cmap = 'RdBu_r'
    if cmass is not None:
        ccmap = 'Blues'
    if kind == 'auto':
        cf = pmass.plot.contourf(ax=ax,levels=levels,x='time',robust=True,zorder=1,cmap=cmap,extend='both')
    elif kind == 'contour':
        if colr is None:
            colr = 'k'
        pmass.plot.contour(ax=ax,levels=levels,x='time',colors=[colr],negative_linestyles='dashed',zorder=1,add_colorbar=False)#,cmap='PuOr')
    if ttle is None:
        ttle = 'Water vapor mass'
    outFile = 'figures/tQ'
    outFile = outFile+appendix+'.pdf'
    ax.set_axisbelow(False)
    fig.subplots_adjust(hspace=0.5)
    #ax = plt.gca()
    ax.set_title(ttle)
    ax.set_xlabel('time [years since eruption]')
    ax.grid(True,zorder=5)
    import datetime
    lim1 = datetime.date(2022,1,15)
    lim2 = datetime.date(2023,11,15)
    ax.set_xlim(lim1,lim2)
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    if fig_out:
        return fig,ax
    else:
        fc.SaveFig(fig,outFile)

    
if args.level is None:
    pure_anom = True
    units = 'mg/m2'
    mls_var = 'anom_hm'
else:
    pure_anom = False
    units = 'ppmv'
    mls_var = 'h2o'


mls = fc.ReadMLS(pure_anom,args.level,adjust_time=False)
if 'time' in mls:
    #anom_hm = mls[mls_var].resample(time='1M',label='left',loffset='14D').mean()
    anom_hm = mls[mls_var]
    if args.level is not None:
        #anom_hm = anom_hm.groupby('time.month') - mls[mls_var+'_clim']
        clim_hm = mls[mls_var+'_clim']
        clim_hm = xr.concat([clim_hm.isel(month=11).assign_coords(month=0),clim_hm,clim_hm.isel(month=0).assign_coords(month=13)],'month')
        clim_hm = clim_hm.assign_coords(month=365/12*clim_hm.month-15)
        clim_hm = clim_hm.interp(month=np.arange(366))
        clim_hm = clim_hm.rename(month='dayofyear')
        anom_hm = anom_hm.groupby('time.dayofyear') - clim_hm
else:
    anom_hm = mls[mls_var] - mls[mls_var+'_clim']
if args.level is None:
    anom_hm = anom_hm*1e3 #mls data is in g/m2 for integrated wv, and ppmv else
else:
    anom_hm = anom_hm*1e6 #mls data is in ppmv
anom_hm.attrs['units'] = units


if 'time' in anom_hm.dims:
    appendix = '_MLS_'
    if args.level is not None:
        appendix = appendix + '{0}hPa'.format(args.level)
    PlotMass(anom_hm,None,appendix,levels=10,ttle='specific humidity at {0}hPa'.format(args.level))
else:
    yr_mnth = [int(a) for a in args.level.split('-')]
    fig,ax = plt.subplots()
    anom_hm.sel(pres=slice(100,1)).plot.contourf(levels=40,x='lat',ax=ax)
    mls_yr = 2021+yr_mnth[0]
    ax.set_title('MLS {0:04d}-{1:02d}'.format(mls_yr,yr_mnth[1]))
    ac.LogPlot(ax[0])
    ac.LogPlot(ax[1])
    ax.set_ylabel('')
    outFile = 'figures/zQ_MLS_{0}.pdf'.format(model_yr)
    fig.savefig(outFile,bbox_inches='tight',transparent=True)
    print(outFile)

