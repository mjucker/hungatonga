#!env python

import xarray as xr
from aostools import climate as ac
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from hungatonga import functions as fc
import calendar
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--qbo',dest='qbo',default=None,help='Only take ensemble members which start in given QBO phase. Either + or - if given.')
parser.add_argument('--years',dest='years',default=None,help='Only show this range of years.')
parser.add_argument('--ylim',dest='ylim',default=None,nargs=2,help='Fix y limits of plot.')
args = parser.parse_args()

option = 'area' #'DU','area'

hemispheres = ['SH'] #['SH','NH']

operation = 'mean' #'mean','extr'

do_spaghetti = False

if args.years is not None:
    years = [int(y) for y in args.years.split(',')]
    yslice = {'time':slice('{:04d}'.format(years[0]),'{:04d}'.format(years[1]))}
if args.ylim is not None:
    ylims = [float(y) for y in args.ylim]

#pert = xr.open_dataset('waccm_pert_ens.nc').TCO
pert = xr.open_dataset('waccm_pert_ens.nc',decode_times=False)
pert,_,_ = fc.CorrectTime(pert)
if args.years is not None:
    pert = pert.sel(yslice)
#ctrl = xr.open_dataset('waccm_ctrl_ens.nc').TCO
ctrl = xr.open_dataset('waccm_ctrl_ens.nc',decode_times=False)
ctrl,_,_ = fc.CorrectTime(ctrl)
if args.years is not None:
    ctrl = ctrl.sel(yslice)
if args.qbo is not None:
    qbo_pos,qbo_neg = fc.CheckQBO(pert,'waccm')
    if args.qbo == '+':
        pert = pert.isel(member=qbo_pos)
        ctrl = ctrl.isel(member=qbo_pos)
    elif args.qbo == '-':
        pert = pert.isel(member=qbo_neg)
        ctrl = ctrl.isel(member=qbo_neg)
pert = pert.TCO
ctrl = ctrl.TCO

if option == 'area':
    oh_pert = {'SH': ac.ComputeOzoneHoleArea(pert)*1e-12, #convert to million km2
               'NH': ac.ComputeOzoneHoleArea(pert,'N')*1e-12}
    oh_ctrl = {'SH' :ac.ComputeOzoneHoleArea(ctrl)*1e-12,
               'NH': ac.ComputeOzoneHoleArea(ctrl,'N')*1e-12}
elif option == 'DU':
    oh_pert = {'SH': ac.GlobalAvgXr(pert.mean('lon'),[-90,-60]),
               'NH': ac.GlobalAvgXr(pert.mean('lon'),[70,90])}
    oh_ctrl = {'SH': ac.GlobalAvgXr(ctrl.mean('lon'),[-90,-60]),
               'NH': ac.GlobalAvgXr(ctrl.mean('lon'),[70,90])}

ylabel = {'area':'area [million km2]',
          'DU'  :'TCO [DU]'
          }


# plot anomalies only
with sns.axes_style('whitegrid'):
    #fig,axs = plt.subplots(ncols=2,figsize=[8,3])
    fig,ax = plt.subplots()
import pandas as pd
max_diff = {}
max_ctrl = {}
months = {'SH':[9,14],'NH':[1,4]}
month_names = {}
for key in hemispheres:
    tf = months[key]
    month_names[key] = [calendar.month_abbr[(t-1)%12+1] for t in tf]
    do_shift = False
    if max(tf) > 12:
        shift = 12-max(tf)
        do_shift = True
    if min(tf) < 0:
        shift = -min(tf)
        do_shift = True
    if do_shift:
        oh_pert[key] = oh_pert[key].shift(time=shift)
        oh_ctrl[key] = oh_ctrl[key].shift(time=shift)
        tf = [t+shift for t in tf]
    filtr = (oh_pert[key]['time.month'] >= months[key][0])*(oh_pert[key]['time.month'] <= months[key][1])
    max_diff[key] = (oh_pert[key]-oh_ctrl[key]).isel(time=filtr)
    if option == 'area':
        if operation == 'mean':
            max_diff[key] = max_diff[key].groupby('time.year').mean()
            op_name = 'mean'
        elif operation == 'extr':
            max_diff[key] = max_diff[key].groupby('time.year').max()
            op_name = 'max'
    elif option == 'DU':
        if operation == 'mean':
            max_diff[key] = max_diff[key].groupby('time.year').mean()
            op_name = 'mean'
        elif operation == 'extr':
            max_diff[key] = max_diff[key].groupby('time.year').min()
            op_name = 'min'
        
colrs = sns.color_palette()
memb = []
years = []
vals = []
hemi = []
for hem in max_diff.keys():
    hem_nme = '{0} {1}-{2}'.format(hem,*month_names[hem])
    for m in range(len(max_diff[hem].member)):
        for y in range(len(max_diff[hem].year)):
            hemi.append(hem_nme)
            memb.append(max_diff[hem].member.values[m])
            years.append(max_diff[hem].year.values[y])
            vals.append(float(max_diff[hem].isel(member=m,year=y).values))
yr1 = years[0]
years = [y-yr1+1 for y in years]
dg = pd.DataFrame(data={'hemi':hemi,'member':memb,'year':years,'TCO':vals})
with sns.axes_style('whitegrid'):
    try: # new sns versions
        sns.barplot(data=dg,x='year',y='TCO',hue='hemi',errorbar=('ci',90),ax=ax)
    except: # old sns vesions
        sns.barplot(data=dg,x='year',y='TCO',hue='hemi',ci=90,ax=ax)
    sns.despine(bottom=True)
if option == 'area':
    ttle = op_name+' ozone hole area [million km2]'
elif option == 'DU':
    ttle = op_name+' polar cap TCO anomaly [DU]'
if args.qbo is not None:
    if args.qbo == '+':
        ttle = ttle+', QBOW'
    elif args.qbo == '-':
        ttle = ttle+', QBOE'
ax.set_title(ttle)
ax.set_ylabel(ylabel[option])
if args.ylim is not None:
    ax.set_ylim(*ylims)
ax.set_xlabel('time [years since eruption]')
outFile = 'figures/tco_anomaly_{0}.pdf'.format(option)
if args.qbo is not None:
    outFile = fc.RenameQBOFile(outFile,args.qbo)
fig.savefig(outFile,transparent=True,bbox_inches='tight')
print(outFile)


if do_spaghetti:
    # Total polar cap ozone

    def MembStats(ds):
        ds_min   = ds.min('member')
        ds_max   = ds.max('member')
        ds_quants= ds.quantile([0.1,0.3,0.7,0.9],'member')
        ds_mean  = ds.mean('member')
        return ds_min,ds_max,ds_mean,ds_quants


    hemis = {'SH':[-90,-63],'NH':[63,90]}

    pyears = np.unique(ctrl.time.dt.year)
    nyears = len(pyears)

    for hemi,sel in hemis.items():
        ctrl_pc = ac.GlobalAvgXr(ctrl.mean(['lon']),sel)
        pert_pc = ac.GlobalAvgXr(pert.mean('lon'),sel)


        fig,axs = plt.subplots(nrows=nyears//2,ncols=2,sharex=True)
        fig.set_figheight(3*nyears//2,4*2)
        ctrl_col = 'k'
        pert_col = 'r'
        for y,year in enumerate(pyears):
            filtr = ctrl_pc['time.year'] == year
            ctrl_min,ctrl_max,ctrl_mean,ctrl_quants = MembStats(ctrl_pc.isel(time=filtr).groupby('time.month').mean())
            ax = axs.flatten()[y]
            p = ctrl_min.plot.line(color=ctrl_col,ax=ax,lw=1)
            ctrl_max.plot.line(color=ctrl_col,ax=ax,lw=1)
            xdata = p[0].get_xdata()
            ax.fill_between(xdata,ctrl_quants.isel(quantile=0),ctrl_quants.isel(quantile=-1),color=ctrl_col,alpha=0.3)
            ax.fill_between(xdata,ctrl_quants.isel(quantile=1),ctrl_quants.isel(quantile=-2),color=ctrl_col,alpha=0.3)
            ctrl_mean.plot.line(color='k',lw=2,ax=ax)
            #
            pert_min,pert_max,pert_mean,pert_quants = MembStats(pert_pc.isel(time=filtr).groupby('time.month').mean())
            p = pert_min.plot.line(color=pert_col,ax=ax)
            pert_max.plot.line(color=pert_col,ax=ax)
            xdata = p[0].get_xdata()
            ax.fill_between(xdata,pert_quants.isel(quantile=0),pert_quants.isel(quantile=-1),color=pert_col,alpha=0.3)
            ax.fill_between(xdata,pert_quants.isel(quantile=1),pert_quants.isel(quantile=-2),color=pert_col,alpha=0.3)
            pert_mean.plot.line(color=pert_col,lw=2,ax=ax)
            ax.set_title('year {0}'.format(y+1))
            sns.despine(ax=ax,top=True,right=True)
            if y < nyears-2:
                ax.set_xlabel('')
            ax.grid()
        ax.set_xticks(ctrl_min.month.values)
        ax.set_xlim(ctrl_min.month.min().values,ctrl_min.month.max().values)
        outFile = 'figures/polar_cap_tco_{0}.pdf'.format(hemi)
        if args.qbo is not None:
            outFile = fc.RenameQBOFile(outFile,args.qbo)
        fig.savefig(outFile,bbox_inches='tight',transparent=True)
        print(outFile)
