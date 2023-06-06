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
clrs = sns.color_palette()
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
args = parser.parse_args()
model = args.model
qmodel = fc.ModName(model)
Q = fc.variables[qmodel]['Q']
U = fc.variables[qmodel]['U']

limits = [1,100]

base = xr.open_dataset(model+'_pert_ens.nc',decode_times=False)
base = xr.merge([base[Q],base[U]])
diff = xr.open_dataset(model+'_delta_ens.nc',decode_times=False)
diff = xr.merge([diff[Q],diff[U]])

diff,_,_ = fc.CorrectTime(diff)
base,_,_ = fc.CorrectTime(base)
dQ = diff[Q]

dimNames = ac.FindCoordNames(dQ)
lev = dimNames['pres']
lat = dimNames['lat']
if dQ[lev][0] > dQ[lev][-1]:
    limits = limits[::-1]

sphum = dQ.sel({lev:slice(*limits)})
mass  = ac.GlobalMass(sphum)/1e9
qbo   = base[U].sel(lat=slice(-5,5)).mean(['lat','lon'])

qbo.plot(x='time',yincrease=False,yscale='log',col='member',col_wrap=5)
outFile = 'figures/{0}_QBO_members.pdf'.format(model)
plt.gcf().savefig(outFile,transparent=True,bbox_inches='tight')
print(outFile)


phase = qbo.sel({lev:50})

pos = phase.isel(time=0) > 0
neg = phase.isel(time=0) < 0
del pos[lev]
del pos['time']
del neg[lev]
del neg['time']

fig,ax = plt.subplots()

mass.mean('member').plot(x='time',ax=ax,color=clrs[0],label='all')
mass.isel(member=pos).mean('member').plot(x='time',ax=ax,color=clrs[1],label='pos')
mass.isel(member=neg).mean('member').plot(x='time',ax=ax,color=clrs[2],label='neg')
ax.legend()

ax2 = ax.twinx()
phase.isel(member=pos).mean('member').plot(x='time',ax=ax2,color=clrs[1],ls='--')
phase.isel(member=neg).mean('member').plot(x='time',ax=ax2,color=clrs[2],ls='--')
ax2.set_title('')
ax2.grid(False)

ax.set_title(model+' total Q by QBO')

outFile = 'figures/{0}_QBO_tQ.pdf'.format(model)
fig.savefig(outFile,transparent=True,bbox_inches='tight')
print(outFile)

fig2,ax2 = plt.subplots()
dphase = phase.isel(member=pos).mean('member') - phase.isel(member=neg).mean('member')
dmass = mass.isel(member=pos).mean('member') - mass.isel(member=neg).mean('member')
dphase.plot(ax=ax2,label='QBO')
dmass.plot(ax=ax2,label='dH2O')
ax2.axhline(0,color='k')
ax2.legend()
ax2.set_ylabel(u'$\Delta$U50 [m/s], $\Delta$M [Tg]')
sns.despine(ax=ax2,bottom=True)
ax2.set_title(model+' QBO, total Q difference WQBO-EQBO')
outFile = 'figures/{0}_dQBO_dtQ.pdf'.format(model)
fig2.savefig(outFile,transparent=True,bbox_inches='tight')
print(outFile)

