#!env python

import xarray as xr
from aostools import climate as ac
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from hungatonga import functions as fc

#pert = xr.open_dataset('waccm_pert_ens.nc').TCO
pert = xr.open_dataset('waccm_pert_ens.nc',decode_times=False)
tmp,_,_ = fc.CorrectTime(pert)
pert = tmp.TCO
#ctrl = xr.open_dataset('waccm_ctrl_ens.nc').TCO
ctrl = xr.open_dataset('waccm_ctrl_ens.nc',decode_times=False)
tmp,_,_ = fc.CorrectTime(ctrl)
ctrl = tmp.TCO

#oh_pert = {'SH': ac.ComputeOzoneHoleArea(pert)*1e-12, #convert to million km2
#           'NH': ac.ComputeOzoneHoleArea(pert,'N')*1e-12}
#oh_ctrl = {'SH' :ac.ComputeOzoneHoleArea(ctrl)*1e-12,
#           'NH': ac.ComputeOzoneHoleArea(ctrl,'N')*1e-12}

oh_pert = {'SH': ac.GlobalAvgXr(pert.mean('lon'),[-90,-60]),
           'NH': ac.GlobalAvgXr(pert.mean('lon'),[70,90])}
oh_ctrl = {'SH': ac.GlobalAvgXr(ctrl.mean('lon'),[-90,-60]),
           'NH': ac.GlobalAvgXr(ctrl.mean('lon'),[70,90])}



# plot anomalies only
with sns.axes_style('whitegrid'):
    #fig,axs = plt.subplots(ncols=2,figsize=[8,3])
    fig,ax = plt.subplots()
import pandas as pd
max_diff = {}
max_ctrl = {}
months = {'SH':[9,9],'NH':[3,3]}
for key in ['NH','SH']:
    #max_diff[key] = oh_pert[key].groupby('time.year').max() - oh_ctrl[key].groupby('time.year').max()
    #max_diff[key] = (oh_pert[key]-oh_ctrl[key]).groupby('time.year').max()
#    filtr = oh_pert[key]['time.month'] == months[key]
    filtr = (oh_pert[key]['time.month'] >= months[key][0])*(oh_pert[key]['time.month'] <= months[key][1])
    max_diff[key] = (oh_pert[key]-oh_ctrl[key]).isel(time=filtr)
    max_diff[key] = max_diff[key].groupby('time.year').mean()
    #max_ctrl[key] = oh_ctrl[key].isel(time=filtr).groupby('time.year').mean()
colrs = sns.color_palette()
memb = []
years = []
vals = []
hemi = []
for hem in max_diff.keys():
    for m in range(len(max_diff[hem].member)):
        for y in range(len(max_diff[hem].year)):
            hemi.append(hem)
            memb.append(max_diff[hem].member.values[m])
            years.append(max_diff[hem].year.values[y])
            vals.append(float(max_diff[hem].isel(member=m,year=y).values))
yr1 = years[0]
years = [y-yr1+1 for y in years]
dg = pd.DataFrame(data={'hemi':hemi,'member':memb,'year':years,'TCO':vals})
with sns.axes_style('whitegrid'):
        sns.barplot(data=dg,x='year',y='TCO',hue='hemi',ci=90,ax=ax)
        sns.despine(bottom=True)
ax.set_title('Polar cap mean TCO anomaly [DU]')
ax.set_ylabel('TCO anomaly [DU]')
outFile = 'figures/tco_anomaly.pdf'
fig.savefig(outFile,transparent=True,bbox_inches='tight')
print(outFile)

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
    fig.savefig(outFile,bbox_inches='tight',transparent=True)
    print(outFile)
