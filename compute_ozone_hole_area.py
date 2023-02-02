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

oh_pert = {'SH': ac.ComputeOzoneHoleArea(pert)*1e-12, #convert to million km2
           'NH': ac.ComputeOzoneHoleArea(pert,'N')*1e-12}
oh_ctrl = {'SH' :ac.ComputeOzoneHoleArea(ctrl)*1e-12,
           'NH': ac.ComputeOzoneHoleArea(ctrl,'N')*1e-12}

oh_2022 = 21.5 # from ozonewatch, September 2022

#fig,ax = plt.subplots()
#for y,year in enumerate(np.unique(oh_pert['time.year'])):
#    ysel = slice('{0:04d}-02-01'.format(year),'{0:04d}-12-31'.format(year))
#    ptmp = oh_pert.sel(time=ysel).mean('member')
#    if y == 0:
#        times = ptmp.time.values
#    ctmp = oh_ctrl.sel(time=ysel).mean('member')
#    ax.plot(times,ctmp.values,color='k',alpha=0.5)
#    ax.plot(times,ptmp.values,color='r',alpha=0.5)

#oh_ctrl.groupby('time.year').max().plot.line(x='year',color='k',alpha=0.5)
#oh_pert.groupby('time.year'),max().plot.line(x='year',color='r',alpha=0.5)

##next idea: seaborn barplots for ctrl and pert, error bars across members
#yr_size = {'ctrl':oh_ctrl.groupby('time.year').max()}
#yr_size['pert'] = oh_pert.groupby('time.year').max()
## plot both ctrl and pert (but note that ctrl is always the same)
#memb = []
#years = []
#vals = []
#cases = []
#for m in range(len(yr_ctrl.member)):
#    for y in range(len(yr_ctrl.year)):
#        for case in yr_size.keys():
#            cases.append(case)
#            years.append(yr_ctrl.year.values[y])
#            memb.append(yr_ctrl.member.values[m])
#            vals.append(float(yr_size[case].isel(member=m,year=y).values))
#yr1 = years[0]
#years = [y-yr1+1 for y in years]
#df = pd.DataFrame(data={'member':memb,'year':years,'case':cases,'size':vals})
#sns.barplot(data=df,x='year',y='size',hue='case')

# plot anomalies only
with sns.axes_style('whitegrid'):
    #fig,axs = plt.subplots(ncols=2,figsize=[8,3])
    fig,ax = plt.subplots()
axs = [ax]
import pandas as pd
max_diff = {}
max_ctrl = {}
months = {'SH':[9,9],'NH':[4,4]}
for key in ['SH']:#,'NH']:
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
dg = pd.DataFrame(data={'hemi':hemi,'member':memb,'year':years,'size':vals})
with sns.axes_style('whitegrid'):
    for h,hem in enumerate(max_diff.keys()):
        filtr = dg['hemi'] == hem
        sns.barplot(data=dg.loc[filtr],x='year',y='size',color=colrs[h],ci=90,ax=axs[h])
        sns.despine(bottom=True)
        #if hem == 'SH':
        #    axs[h].text(1,oh_2022-max_ctrl[hem].mean('member').isel(year=0).values,'x',color='r')
        #axs[h].set_title(hem)
        #if h == 0:
        axs[h].set_title('Antarctic ozone hole area anomaly [1e6 km2]')
        #else:
        axs[h].set_ylabel('Ozone hole area [1e6 km2]')
outFile = 'figures/ozone_hole_anomaly.pdf'
fig.savefig(outFile,transparent=True,bbox_inches='tight')
print(outFile)
