#!env python

import xarray as xr
from aostools import climate as ac
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from hungatonga import functions as fc
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
parser.add_argument('--qbo',dest='qbo',default=None,help='Only take ensemble members which start in given QBO phase. Either + or - if given.')
parser.add_argument('--qmodel',dest='qbo_model',default=None,help='Use this model to assess QBO phase.')
args = parser.parse_args()
model = args.model
qnme = {'+':'QBOW','-':'QBOE'}

sns.set_context('notebook')
sns.set_style('whitegrid')

def ReadFile(file):
    ds = xr.open_dataset(file,decode_times=False)
    ds,units,cal = fc.CorrectTime(ds,decode=False)
    attrs = ds.time.attrs
    if cal.lower() == 'noleap':
        daysperyear = 365
    else:
        daysperyear = 360
    ds = ds.assign_coords(time=(ds.time/daysperyear))
    ds.time.attrs['units'] = 'years'
    #ds.time.attrs['calendar'] = attrs['calendar']
    return xr.decode_cf(ds)

qmodel = fc.ModName(model)
Q = fc.variables[qmodel]['Q']
pert_file = model+'_pert_ens.nc'
ctrl_file = model+'_ctrl_ens.nc'
pert = ReadFile(pert_file)
ctrl = ReadFile(ctrl_file)
lev = ac.FindCoordNames(pert)['pres']
if args.qbo is not None:
    if args.qbo_model is None:
        qbo_pos,qbo_neg = fc.CheckQBO(pert,model)
    else:
        pertq = xr.open_dataset(args.qbo_model+'_pert_ens.nc',decode_times=False)
        qbo_pos,qbo_neg = fc.CheckQBO(pertq,args.qbo_model)
        qbo_pos = qbo_pos.assign_coords({'member':pert.member})
        qbo_neg = qbo_neg.assign_coords({'member':pert.member})
    if args.qbo == '+':
        pert = pert.isel(member=qbo_pos)
        ctrl = ctrl.isel(member=qbo_pos)
    elif args.qbo == '-':
        pert = pert.isel(member=qbo_neg)
        ctrl = ctrl.isel(member=qbo_neg)
if model.lower() == 'waccm':
    levslice = {lev:slice(1,100)}
else:
    levslice = {lev:slice(100,1)}
pert = pert.sel(levslice)
ctrl = ctrl.sel(levslice)
if not 'CLDICE' in pert.data_vars:
    pert['CLDICE'] = 0*pert[Q]
    ctrl['CLDICE'] = 0*ctrl[Q]
pert = xr.merge([pert[Q],pert.CLDICE])
ctrl = xr.merge([ctrl[Q],ctrl.CLDICE])



do_daily = False
pert_file_d = 'waccm_daily_pert_ens.nc'
if os.path.isfile(pert_file_d):
    do_daily = True
    pert_d = ReadFile(pert_file_d).sel(levslice)
    pert_d = xr.merge([pert_d.CLDICE,pert_d.Q])

def MassPerDeg(ds):
    from aostools.constants import a0,g,coslat
    mass_y = coslat(ds.lat)*ds
    mass_x = np.deg2rad(mass_y.integrate('lon'))
    mass_p = mass_x.integrate(lev)*100
    mass = a0**2/g*mass_p
    if ds[lev][0] > ds[lev][-1]:
        mass = -mass
    return mass


def PlotMass(mass,tco=None,appendix='',levels=20,colr=None,fig_out=False,kind='auto',fig=None,ax=None):
    #fig,ax = plt.subplots()
    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = fig.add_axes([0.1,0.1,0.74,0.8])
    cmass = None
    if isinstance(mass,xr.Dataset):
        if Q in mass.data_vars:
            qmass = mass[Q]
        if 'CLDICE' in mass.data_vars:
            cmass = mass.CLDICE
    else:
        qmass = mass
    if 'member' in qmass.coords:
        pval = ac.StatTest(qmass,0,'T','member',parallel=True)
        pmass = qmass.mean('member').where(pval<0.1)
        smass = qmass.mean('member')
    else:
        pmass = qmass
    if cmass is not None and 'member' in cmass.coords:
        cval = ac.StatTest(cmass,0,'T','member',parallel=True)
        cmass = cmass.mean('member').where(cval<0.1)
    if qmass.name == Q:
        cmap = 'Reds'
    else:
        cmap = 'RdBu_r'
    if cmass is not None:
        ccmap = 'Blues'
    if kind == 'auto':
        cf = pmass.plot.contourf(ax=ax,levels=levels,x='time',robust=True,zorder=1,cmap=cmap,add_colorbar=False)#,cbar_kwargs={'label':'Q [mg/m2]','shrink':0.8},cmap=cmap)#,cmap='PuOr')
        position = fig.add_axes([0.85,0.12,0.02,0.35])
        fig.colorbar(cf,ax=ax,cax=position,label='Q [mg/m2]')
        if cmass is not None:
            cf = cmass.plot.contourf(ax=ax,levels=np.linspace(0,50,11),x='time',robust=True,zorder=2,add_colorbar=False,cmap=ccmap)#,cbar_kwargs={'format':'%4.1f','label':'CLDICE [mg/m2]','shrink':0.5,'anchor':(0,0.75)})
            position = fig.add_axes([0.85,0.5,0.02,0.35])
            fig.colorbar(cf,ax=ax,cax=position,label='CLDICE [mg/m2]')
    elif kind == 'contour':
        if colr is None:
            colr = 'k'
        pmass.plot.contour(ax=ax,levels=levels,x='time',colors=colr,zorder=1,add_colorbar=False)#,cmap='PuOr')
    ttle = 'Water vapor mass'
    outFile = 'figures/tQ'
    if tco is not None:
        pval = ac.StatTest(tco,0,'T','member',parallel=True)
        tco.mean('member').plot.contour(ax=ax,levels=[-10,-5],colors='k',x='time',zorder=2,linestyles='dashed')
        tco.mean('member').where(pval<0.1).plot.contour(ax=ax,levels=[-10,-5],colors='k',x='time',zorder=3,linestyles='solid')
        ttle = ttle+', TCO anomalies'
        outFile = outFile+'_dTCO'
    if args.qbo is not None:
        ttle = ttle+', '+qnme[args.qbo]
    outFile = outFile+appendix+'.pdf'
    if args.qbo is not None:
        outFile = fc.RenameQBOFile(outFile,args.qbo)
    ax.set_axisbelow(False)
    fig.subplots_adjust(hspace=0.5)
    #ax = plt.gca()
    ax.set_title(ttle)
    ax.set_xlabel('time [years since eruption]')
    ax.grid(True,zorder=5)
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    if fig_out:
        return fig,ax
    else:
        fc.SaveFig(fig,outFile)



delta = fc.Mass(pert-ctrl)*1e3
#delta = delta.mean('member')
if isinstance(delta,xr.Dataset):
    for var in delta.data_vars:
        delta[var].attrs['units'] = 'mg/m2'
else:
    delta.attrs['units'] = 'g/m2'
if do_daily:
    delta_d = fc.Mass(pert_d)*1e3
    delta_d = delta_d.mean('member')
    delta_d.attrs['units'] = 'g/m2'

levs = [2,4,6,8,10]

mls = fc.ReadMLS(True).sel(time=slice('0001-01-01',None))
anom_hm = mls.anom_hm.resample(time='1M',label='left',loffset='14D').mean()
anom_hm = fc.ConvertTime2Years(anom_hm)

## Ozone
if model.lower() == 'waccm':
    dTCO = ReadFile(pert_file).TCO - ReadFile(ctrl_file).TCO
    dTCO = dTCO.mean('lon')
else:
    dTCO = None
if args.qbo is not None and dTCO is not None:
    if args.qbo == '+':
        dTCO = dTCO.isel(member=qbo_pos)
    elif args.qbo == '-':
        dTCO = dTCO.isel(member=qbo_neg)

fig,ax = PlotMass(delta.mean('lon')*1e3,dTCO,'_'+model.upper(),fig_out = True)
PlotMass(anom_hm,None,'_{0}_MLS'.format(model.upper()),kind='contour',colr='cyan',levels=[.25,.5,.75,1],fig=fig,ax=ax)

#delta = delta.mean('member')

# 2D early daily evolution
if do_daily:
    mls = fc.ReadMLSMap()
    mls = ac.CloseGlobe(mls)
    mls = mls - mls.sel(time=slice(None,'0001-01-14')).mean('time')

    cjan = fc.Mass(ctrl[Q].isel(time=0).mean('member'))*1e3
    inds = [14,21,28,35]
    nrows = 2
    fig,axs,transf = ac.Projection('Robinson',ncols=int(len(inds)/nrows),nrows=nrows,kw_args={'central_longitude':180})
    fig.set_figwidth(6*len(inds)/nrows)

    for a,ax in enumerate(axs.flatten()):
        mls.isel(time=inds[a]).plot.contour(ax=ax,levels=levs,cmap='viridis',**transf)
        ttle = ax.get_title().split(' ')[2]
        ttle = ttle.replace('0001','2022')
        cf = (delta_d[Q]-cjan).isel(time=inds[a]-inds[0]+1).plot(levels=11,vmin=0,vmax=3,extend='max',ax=ax,cmap='Blues',add_colorbar=False,**transf)
        ax.set_title(ttle)
        ax.gridlines()
        ax.coastlines()
    ac.AddColorbar(fig,axs,cf,shrink=0.6,cbar_args={'label':'Q [g/m2]'})
    ac.AddPanelLabels(axs,'upper left',ypos=1.1)
    fc.SaveFig(fig,'figures/cloud_days.pdf')

    # Ice cloud formation in WACCM
    nrows = 2
    fig,axs,transf = ac.Projection('Robinson',ncols=2,nrows=nrows,kw_args={'central_longitude':180})
    fig.set_figwidth(6*2)

    for a,ax in enumerate(axs.flatten()):
        cf = (delta_d.CLDICE).isel(time=a).plot(levels=13,vmin=0,vmax=1.2,extend='max',ax=ax,cmap='gray',add_colorbar=False,**transf)
        ax.set_title('day {0}'.format(a))
        #ax.gridlines()
        ax.coastlines(color='gray')
    ac.AddColorbar(fig,axs,cf,shrink=0.6,cbar_args={'label':'CLDICE [g/m2]'})
    ac.AddPanelLabels(axs,'upper left',ypos=1.1)
    fc.SaveFig(fig,'figures/ice_clouds.pdf')
