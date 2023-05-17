#!env python
import xarray as xr
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
from matplotlib import pyplot as plt
import numpy as np
sns.set_context('paper')
sns.set_style('whitegrid')
sns.set_palette('tab10')
colrs = sns.color_palette()
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
args = parser.parse_args()
model = args.model
qmodel = fc.ModName(model)
Q = fc.variables[qmodel]['Q']

limits = [1,100]

same_fig = False

pert = xr.open_dataset(model+'_pert_ens.nc',decode_times=False)#.Q
if 'CLDICE' not in pert.data_vars:
    pert['CLDICE'] = pert[Q]*0
pert = xr.merge([pert[Q],pert.CLDICE])
ctrl = xr.open_dataset(model+'_ctrl_ens.nc',decode_times=False)#.Q
if 'CLDICE' not in ctrl.data_vars:
    ctrl['CLDICE'] = pert[Q]*0
ctrl = xr.merge([ctrl[Q],ctrl.CLDICE])
mls = fc.ReadMLS(pure_anom=True)
mls  = mls.sel(time=slice('0001-01-01',None))
tq_mls = mls.anom

dTS = pert - ctrl
#dTS,_,_ = fc.CorrectTime(dTS.to_dataset())
dTS,_,_ = fc.CorrectTime(dTS)
dCI = dTS.CLDICE
dTS = dTS[Q]

dimNames = ac.FindCoordNames(dTS)
lev = dimNames['pres']
lat = dimNames['lat']
if dTS[lev][0] > dTS[lev][-1]:
    limits = limits[::-1]

sphum = dTS.sel({lev:slice(*limits)})
cldiceS= dCI.sel({lat:slice(-90,-60),lev:slice(*limits)})
cldiceN= dCI.sel({lat:slice(60,90),lev:slice(*limits)})
tq = ac.GlobalMass(sphum)
tcS = ac.GlobalMass(cldiceS)
tcN = ac.GlobalMass(cldiceN)
if same_fig:
    fig,axs = plt.subplots(nrows=2,sharex=True)
    labl1 = 0.6
    labl2 = 0.43
else:
    labl1 = 0.7
    labl2 = 0.63
    fig,ax = plt.subplots()
    fig2,ax2 = plt.subplots()
    axs = [ax,ax2]
fc.GlobalMeanPlot(tq*1e-9,name=None,fig=fig,ax=axs[0],color=colrs[3],label=model.upper())
tq_mls_mn = tq_mls.resample(time='1M',label='left',loffset='16D').mean()
p = tq_mls_mn.plot(ax=axs[0],label='MLS',color=colrs[0],lw=2,ls='-')
axs[0].legend()
axs[0].fill_between(p[0].get_xdata(),tq_mls_mn-4,tq_mls_mn+4,color=p[0].get_color(),alpha=0.3)
axs[0].set_title('Total stratospheric water mass [Tg]')
# add total cloud ice
fc.GlobalMeanPlot(tcS*1e-9,name=None,fig=fig,ax=axs[1],color=colrs[2],label='60-90S')
fc.GlobalMeanPlot(tcN*1e-9,name=None,fig=fig,ax=axs[1],color=colrs[4],label='60-90N',ttle='Polar stratospheric cloud ice mass [Tg]')
axs[1].legend()#['WACCM','MLS'])
## exponential fit
eps = 16/365
# first, only the first few months
nmonths = 6
wtime = np.linspace(eps,nmonths/12-eps,nmonths)
tqm = 1e-9*tq.mean('member').isel(time=slice(0,nmonths))
wlifetime,whalftime,ppm,ppy,amp = fc.ExpFit(wtime,tqm)
print(model.upper()+' first {0} months rate = {1:.2f}% per month'.format(nmonths,ppm))
# then, for the longer term
fit_years = 5
wtime = np.linspace(eps,fit_years-eps,12*fit_years)
tqm = 1e-9*tq.mean('member').isel(time=slice(0,12*fit_years))
wlifetime,whalftime,ppm,ppy,amp = fc.ExpFit(wtime,tqm)
print(model.upper()+' rate [percent per month] = {0:6.2f}'.format(ppm))
axs[0].text(0.8,labl1,r'$\tau_{0} = ${1:3.1f} yrs'.format('{1/2}',whalftime),transform=axs[0].transAxes,color='r',backgroundcolor='white')
axs[0].plot(tqm.time,np.exp(amp)*np.exp(-wtime/wlifetime),color='r',lw=1,ls=':')
nt = len(tq_mls_mn.time)
mtime = np.linspace(eps,nt/12-eps,nt)
mlifetime,mhalftime,ppm,ppy,amp = fc.ExpFit(mtime,tq_mls_mn)
print('MLS rate [percent per month] = {0:6.2f}'.format(ppm))
axs[0].text(0.8,labl2,r'$\tau_{0} = ${1:3.1f} yrs'.format('{1/2}',mhalftime),transform=axs[0].transAxes,color=colrs[0],backgroundcolor='white')
axs[0].plot(tq_mls_mn.time,np.exp(amp)*np.exp(-mtime/mlifetime),color=colrs[0],lw=1,ls=':')
sns.despine(bottom=True,left=True)
#plt.grid()
axs[0].grid()
if same_fig:
    axs[0].set_xlabel('')
else:
    axs[0].set_xlabel('time [years since eruption]')
axs[1].set_xlabel('time [years since eruption]')
axs[0].set_ylabel('mass anomaly [Tg]')
axs[1].set_ylabel('mass anomaly [Tg]')
#ax.set_title('Total stratospheric water mass [Tg]')
try:
    tdat = axs[0].get_children()[0].get_xdata()
    axs[0].set_xlim(tdat[0],tdat[-1])
except:
    print('OBS: could not adjust xlims')
fc.SaveFig(fig,'figures/{0}_tQ.pdf'.format(model))
if not same_fig:
    try:
        axs[1].set_xlim(tdat[0],tdat[-1])
    except:
        print('OBS: could not adjust xlims')
    fc.SaveFig(fig2,'figures/{0}_tPSC.pdf'.format(model))


