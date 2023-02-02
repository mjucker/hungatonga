#!env python
import xarray as xr
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator,ScalarFormatter
from shapely import geometry
from cartopy import crs as ccrs
import numpy as np
import calendar
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',help='Choose model/run to plot.')
args = parser.parse_args()
model = args.model
sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()

fields = {'ZM'  : ['Q','T','U','OMEGA'],
          'Maps': ['OLR','P','SLP','TS']
          }
trans_fields = {'TS':'TREFHT'}


sels = {'Q' : {'lev':slice(1,100)},
        'OMEGA':{'lev':slice(100,800)}
}

lins = ['OMEGA']

scales = {'TS': 1, 'OLR': 1, 'Q': 1e6, 'U': 1, 'T': 1, 'OMEGA': 1e2}

if 'waccm' in model:
    scales['P'] = 1e3*86400
    scales['SLP'] = 0.01
else:
    scales['P'] = 86400
    scales['SLP'] = 1

months = [0,1,2]

loncents = {
    'P' : 155,
    'OLR':155,
    'TS': 155,
    'SLP':155,
}

vmins = {
#    'Q' : 0,
    }
vmaxs = {
#    'Q' : 3,
#    'T' : 0.8,
#    'U' : 2.5,
    'TS': 1.5,
    'P' : 2.0,
    'SLP': 3,
    'OLR': 10,
    }

levs = {
    'Q' : np.linspace(0,3,16),
    'T' : [l for l in np.linspace(-1.0,1.0,11) if l != 0],
    'U' : [l for l in np.linspace(-2.4,2.4,13) if l != 0],
    'OMEGA':np.linspace(-0.3,0.3,16),
}

cmaps = {
    'Q' : 'Reds',
    'U' : 'RdYlBu_r',
    'T' : 'RdBu_r',
    'OMEGA':'BrBG_r',
    'TS': 'RdBu_r',
    'P' : 'PuOr',
    'OLR':'PiYG',
    'SLP':'BrBG_r',
    }

labls = {
    'Q' : 'Q [ppm]',
    'U' : 'U [m/s]',
    'T' : 'T [K]',
    'OMEGA': 'w [hPa/s]',
    'TS': 'T2m [K]',
    'P' : 'P [mm/day]',
    'OLR':'OLR [W/m2]',
    'SLP':'SLP [hPa]',
    }

pert = xr.open_dataset('{0}_pert_ens.nc'.format(model),decode_times=False)
ctrl = xr.open_dataset('{0}_ctrl_ens.nc'.format(model),decode_times=False)


dTS = pert - ctrl
dTS,_,_ = fc.CorrectTime(dTS)
dTS = dTS.isel(time=slice(None,max(months)+1))

#psip = xr.open_dataset('{0}_psis_pert_ens.nc'.format(model)).psis
#psic = xr.open_dataset('{0}_psis_ctrl_ens.nc'.format(model)).psis
#psis = (psip - psic).isel(time=slice(None,max(months)+1))

#add_contours = {'T':{'data': psis, 'levels':[-2e9,-1e9]}}

#ep_pert = xr.open_dataset('{0}_ep_pert_ens.nc'.format(model))
#ep_ctrl = xr.open_dataset('{0}_ep_ctrl_ens.nc'.format(model))
#ep_diff = (ep_pert-ep_ctrl).sel(lat=slice(-80,80),lev=slice(1,800)).isel(time=slice(None,max(months)+1)).isel(lat=slice(None,None,2),lev=slice(None,None,2))
#ep_add = ['U']

hatches = ['//']
FigArgs = {'levels':20,'x':'lat','yincrease':False,'zorder':1}
# First, lat-pres plots
nvars = len(fields['ZM'])
nmonths = len(months)
fig,axs = plt.subplots(figsize=[6*nmonths,3*0.8*nvars],ncols=nmonths,nrows=nvars,sharey=False,sharex=True)
for f,field in enumerate(fields['ZM']):
    var = fc.variables[model][field]
    if var in trans_fields.keys():
        var = trans_fields[var]
    da = (dTS[var]*scales[field]).mean('lon')
    if field in sels.keys():
        da = da.sel(sels[field])
    tmpArgs = {}
    for keys,vals in FigArgs.items():
        tmpArgs[keys] = vals
    if field in lins:
        tmpArgs['yscale'] = 'linear'
    else:
        tmpArgs['yscale'] = 'log'
    #if field in ep_add:
    #    tmpArgs['yincrease'] = False
    tmpArgs['levels'] = levs[field]
    if field in vmins.keys():
        tmpArgs['vmin'] = vmins[field]
    if field in vmaxs.keys():
        tmpArgs['vmax'] = vmaxs[field]
    for m,month in enumerate(months):
        pval = ac.StatTest(da.isel(time=month),0,'T','member',parallel=True)
        dm = da.isel(time=month).mean('member')
        #if m == nmonths-1:
        #    addbar = True
        #    tmpArgs['cbar_kwargs'] = {'label':labls[field]}
        #else:
        addbar = False
        cf=dm.plot.contourf(ax=axs[f][m],cmap=cmaps[field],add_colorbar=addbar,**tmpArgs)
        dm.where(pval<0.1).plot.contourf(ax=axs[f][m],colors='none',hatches=hatches,add_colorbar=False,**FigArgs)
        #if field in ep_add:
        #    pval1 = ac.StatTest(ep_diff.isel(time=month).ep1,0,'T','member',parallel=True)
        #    pval2 = ac.StatTest(ep_diff.isel(time=month).ep2,0,'T','member',parallel=True)
        #    pval = pval1 + pval2
        #    ep1tmp = ep_diff.isel(time=month).ep1.mean('member').where(pval<0.1)
        #    ep2tmp = ep_diff.isel(time=month).ep2.mean('member').where(pval<0.1)
         #   if field in sels.keys():
         #       ep1tmp = ep1tmp.sel(sels[field])
         #       ep2tmp = ep2tmp.sel(sels[field])
         #   ac.PlotEPfluxArrows(ep1tmp.lat,ep1tmp.lev,ep1tmp,ep2tmp,fig,axs[f][m],yscale=tmpArgs['yscale'])
        axs[f][m].set_axisbelow(False)
        axs[f][m].grid(True,zorder=2)
        axs[f][m].set_title('{0}, month {1}'.format(field,month+1))
        if m == 0:
            axs[f][m].set_ylabel('pressure [hPa]')
        else:
            axs[f][m].set_ylabel('')
        if f < nvars-1:
            axs[f][m].set_xlabel('')
        else:
            axs[f][m].set_xlabel('latitude')
        axs[f][m].yaxis.set_major_formatter(ScalarFormatter())
        if field not in sels.keys():
            axs[f][m].set_ylim(1000,1)
    ac.AddColorbar(fig,axs[f],cf,cbar_args={'label':labls[field]})
ac.AddPanelLabels(axs,'upper left',ypos=1.1)
outFile = 'figures/{0}_zm_months.pdf'.format(model)
fig.savefig(outFile,transparent=True,bbox_inches='tight')
print(outFile)

hatches = ['///']
nvars = len(fields['Maps'])
fig = plt.figure(figsize=[6*nmonths,3.5*nvars])
transf = {'transform':ccrs.PlateCarree()}

axs = []
p = 0
for f,field in enumerate(fields['Maps']):
    var = fc.variables[model][field]
    if var in trans_fields.keys():
        var = trans_fields[var]
    da = dTS[var]*scales[field]
    axsr = []
    for m,month in enumerate(months):
        p += 1
        ax = fig.add_subplot(nvars,nmonths,p,projection=ccrs.Robinson(central_longitude=loncents[field]))
        axs.append(ax)
        axsr.append(ax)
        seas = calendar.month_abbr[month+1]
        dm = da.isel(time=month)
        pval = ac.StatTest(dm,0,'T','member',parallel=True)
        dmm = dm.mean('member')
        #if m == nmonths-1:
        #    tmpArgs = {'cbar_kwargs':{'label':labls[field],'shrink':0.85}}
        #else:
        tmpArgs = {'add_colorbar':False}
        tmpArgs['transform'] = transf['transform']
        dmm.where(pval<0.1).plot.contourf(ax=ax,colors='none',hatches=hatches,add_colorbar=False,zorder=2,**transf)
        cf=dmm.plot(ax=ax,vmax=vmaxs[field],cmap=cmaps[field],zorder=1,**tmpArgs)
        ax.gridlines()
        ax.coastlines()
        ax.set_title(field+', month {0}'.format(month+1))
    ac.AddColorbar(fig,axsr,cf,shrink=0.8,cbar_args={'label':labls[field]})
ac.AddPanelLabels(axs,'upper left',ypos=1.1) 
outFile = 'figures/{0}_maps_months.pdf'.format(model)
fig.savefig(outFile,bbox_inches='tight',transparent=True)
print(outFile)
