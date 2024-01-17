#!env python
import xarray as xr
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
from matplotlib import pyplot as plt
from datetime import timedelta
import numpy as np
sns.set_context('notebook')
sns.set_style('whitegrid')
sns.set_palette('tab10')
colrs = sns.color_palette()
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
parser.add_argument('--qbo',dest='qbo',default=None,help='Only take ensemble members which start in given QBO phase. Either + or - (for one phase only) or a (all phases plus total) if given.')
parser.add_argument('--qmodel',dest='qbo_model',default=None,help='Use this model to assess QBO phase.')
args = parser.parse_args()
model = args.model
qmodel = fc.ModName(model)
Q = fc.variables[qmodel]['Q']

limits = [1,100]

tq_limits = None

PSC_lims = [-3,3]

same_fig = False

do_exp = True

qnme = {'+':'W','-':'E'}

pert = xr.open_dataset(model+'_pert_ens.nc',decode_times=False)#.Q
if 'CLDICE' not in pert.data_vars:
    pert['CLDICE'] = pert[Q]*0
ctrl = xr.open_dataset(model+'_ctrl_ens.nc',decode_times=False)#.Q
if 'CLDICE' not in ctrl.data_vars:
    ctrl['CLDICE'] = pert[Q]*0
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
pert = xr.merge([pert[Q],pert.CLDICE])
ctrl = xr.merge([ctrl[Q],ctrl.CLDICE])
mls = fc.ReadMLS(pure_anom=True)
mls  = mls.sel(time=slice('0001-01-01',None))
tq_mls = mls.anom
tq_mls = tq_mls.assign_coords({'time':tq_mls.time+timedelta(days=-365)})

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
    labl1 = 0.5
    labl2 = 0.63
    fig,ax = plt.subplots()
    fig2,ax2 = plt.subplots()
    axs = [ax,ax2]
    if args.qbo == 'a':
        fig3,ax3 = plt.subplots()
        axs.append(ax3)

if args.qbo == 'a':
    fc.GlobalMeanPlot(tq.isel(member=qbo_neg)*1e-9,name=None,fig=fig,ax=axs[0],color=colrs[4],label=model.upper()+' QBO'+qnme['-'],ls='--',fill=False)
    fc.GlobalMeanPlot(tq.isel(member=qbo_pos)*1e-9,name=None,fig=fig,ax=axs[0],color=colrs[2],label=model.upper()+' QBO'+qnme['+'],ls='-.',fill=False)
fc.GlobalMeanPlot(tq*1e-9,name=None,fig=fig,ax=axs[0],color=colrs[3],label=model.upper())
tq_mls_mn = tq_mls.resample(time='1M',label='left',loffset='16D').mean()
p = tq_mls_mn.plot(ax=axs[0],label='MLS',color='k',lw=2,ls='-')
axs[0].legend()
axs[0].fill_between(p[0].get_xdata(),tq_mls_mn-4,tq_mls_mn+4,color=p[0].get_color(),alpha=0.3)
ttle = 'Total stratospheric water mass [Tg]'
if args.qbo is not None and args.qbo != 'a':
    ttle = ttle+', QBO'+qnme[args.qbo]
axs[0].set_title(ttle)

#################3
# add total cloud ice
def PlotCloudIce(tcS,tcN,ax,qbo_str=None):
    fc.GlobalMeanPlot(tcS*1e-9,name=None,fig=fig,ax=ax,color=colrs[2],label='60-90S')
    ttle = 'Polar stratospheric cloud ice mass [Tg]'
    if qbo_str is not None:
        ttle = ttle+', QBO'+qnme[qbo_str]
    fc.GlobalMeanPlot(tcN*1e-9,name=None,fig=fig,ax=ax,color=colrs[4],label='60-90N',ttle=ttle)
    ax.legend()#['WACCM','MLS'])

if args.qbo == 'a':
    PlotCloudIce(tcS.isel(member=qbo_neg),tcN.isel(member=qbo_neg),axs[1],'-')
    PlotCloudIce(tcS.isel(member=qbo_pos),tcN.isel(member=qbo_pos),axs[2],'+')
else:
    PlotCloudIce(tcS,tcN,axs[1],args.qbo)
###############
## exponential fit
if do_exp:
    eps = 16/365
    # first, only the first few months
    nmonths = 6
    wtime = np.linspace(eps,nmonths/12-eps,nmonths)
    tqm = 1e-9*tq.mean('member').isel(time=slice(0,nmonths))
    wlifetime,whalftime,ppm,ppy,amp = fc.ExpFit(wtime,tqm)
    print(model.upper()+' first {0} months rate = {1:.2f}% per month'.format(nmonths,ppm))
    # then, for the longer term
    fit_years = 5
    start_fit = 1
    wtime = np.linspace(start_fit+eps,start_fit+fit_years-eps,12*fit_years)
    tqm = 1e-9*tq.mean('member').isel(time=slice(start_fit*12,12*(start_fit+fit_years)))
    wlifetime,whalftime,ppm,ppy,amp = fc.ExpFit(wtime,tqm)
    print(model.upper()+' rate [percent per month] = {0:6.2f}'.format(ppm))
    print(model.upper()+' e-folding time = {0:6.2f} years'.format(wlifetime))
    #axs[0].text(0.8,labl1,r'$\tau_{0} = ${1:3.1f} yrs'.format('{1/2}',whalftime),transform=axs[0].transAxes,color='r',backgroundcolor='white')
    axs[0].text(0.8,labl1,r'$\tau = ${0:3.1f} yrs'.format(wlifetime),transform=axs[0].transAxes,color=colrs[3],backgroundcolor='white')
    axs[0].plot(tqm.time,np.exp(amp)*np.exp(-wtime/wlifetime),color=colrs[3],lw=1.5,ls=':')
    # MLS residence time
    #nt = len(tq_mls_mn.time)
    #mtime = np.linspace(eps,nt/12-eps,nt)
    #mlifetime,mhalftime,ppm,ppy,amp = fc.ExpFit(mtime,tq_mls_mn)
    #print('MLS rate [percent per month] = {0:6.2f}'.format(ppm))
    #axs[0].text(0.8,labl2,r'$\tau_{0} = ${1:3.1f} yrs'.format('{1/2}',mhalftime),transform=axs[0].transAxes,color=colrs[0],backgroundcolor='white')
    #axs[0].plot(tq_mls_mn.time,np.exp(amp)*np.exp(-mtime/mlifetime),color=colrs[0],lw=1,ls=':')

#######
# figure aesthetics
sns.despine(bottom=True,left=True)
#plt.grid()
axs[0].grid()
if same_fig:
    axs[0].set_xlabel('')
else:
    axs[0].set_xlabel('time [years since eruption]')
for ax in axs[1:]:
    ax.set_xlabel('time [years since eruption]')
for ax in axs:
    ax.set_ylabel('mass anomaly [Tg]')
    ax.set_ylabel('mass anomaly [Tg]')
#ax.set_title('Total stratospheric water mass [Tg]')
try:
    tdat = axs[0].get_children()[0].get_xdata()
    axs[0].set_xlim(tdat[0],tdat[-1])
except:
    print('OBS: could not adjust xlims')
outFile = 'figures/{0}_tQ.pdf'.format(model)
if args.qbo is not None and args.qbo != 'a':
    outFile = fc.RenameQBOFile(outFile,args.qbo)
fc.SaveFig(fig,outFile)
if not same_fig:
    try:
        for ax in axs[1:]:
            ax.set_xlim(tdat[0],tdat[-1])
    except:
        print('OBS: could not adjust xlims')
    if PSC_lims is not None:
        for ax in axs[1:]:
            ax.set_ylim(*PSC_lims)
    outFile = 'figures/{0}_tPSC.pdf'.format(model)
    if args.qbo is None:
        fc.SaveFig(fig2,outFile)
    else:
        if args.qbo == 'a':
            for fig,qbo in zip([fig2,fig3],['-','+']):
                outFileQ = fc.RenameQBOFile(outFile,qbo)
                fc.SaveFig(fig,outFileQ)
        else:
            outFile = fc.RenameQBOFile(outFile,args.qbo)
            fc.SaveFig(fig2,outFile)


