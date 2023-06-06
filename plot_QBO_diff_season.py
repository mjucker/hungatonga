import xarray as xr
import numpy as np
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
from matplotlib.ticker import ScalarFormatter,FormatStrFormatter
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
parser.add_argument('--v1',dest='var1',default='Q',help='Choose the variable to be plotted as shading.')
parser.add_argument('--v2',dest='var2',default='psis',help='Choose the variable to be plotted as contours.')
parser.add_argument('-l',dest='levs',default='1,100',help='Pressure levels to plot.')
parser.add_argument('-D',dest='diff',default='delta',help='Compute the difference between QBO phases of - anomalies [diff,default], - perturbation sims [pert] - control sims [ctrl].')
args = parser.parse_args()
model = args.model
qmodel = fc.ModName(model)
sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()

show = args.diff

plevs = [float(i) for i in args.levs.split(',')]

cutoff_year = 3

vlevels = {
    'O3'  : [l for l in np.linspace(-1e-7,1e-7,21) if l !=0],
    'Q'   : [l for l in np.linspace(-1e-6,1e-6,21) if l !=0],
    'psis': [l for l in np.linspace(-1e11,1e11,21) if l !=0],
    'psi' : [l for l in np.linspace(-1e10,1e10,21) if l !=0],
    }

delta = xr.open_dataset('{0}_season_{1}_ens.nc'.format(model,show),decode_times=False)

if 'psi' in args.var1 or 'psi' in args.var2:
    if 'psi' in args.var1:
        psivar = args.var1
    if 'psi' in args.var2:
        psivar = args.var2
    tmp = xr.open_dataset('{0}_psis_season_{1}_ens.nc'.format(model,show),decode_times=False)
    tmp = tmp[psivar]
    delta = xr.merge([delta,tmp])
delta,_,_ = fc.CorrectTime(delta)
delta = delta.sel(time=slice(None,'{:04d}-11'.format(cutoff_year+1)))
# qbo
ctrl = xr.open_dataset('{0}_ctrl_ens.nc'.format(model),decode_times=False)
qbo_pos,qbo_neg = fc.CheckQBO(ctrl,qmodel)
# zonal means
delta = delta.mean('lon')

# difference of differences depending on QBO phase
dmerge = []
for var in [args.var1,args.var2]:
    pval = ac.StatTest(delta[var].isel(member=qbo_pos),delta[var].isel(member=qbo_neg),'KS','member',parallel=True)
    pval.name = 'pval_'+var
    dmerge.append(pval)
delta = delta.isel(member=qbo_pos).mean('member') - delta.isel(member=qbo_neg).mean('member')
dmerge.insert(0,delta)
delta = xr.merge(dmerge)


#seasonal timeseries, all years
lev = ac.FindCoordNames(delta)['pres']
levs = slice(*plevs)
dtmp = delta.sel({lev:levs})
fg = dtmp[args.var1].where(dtmp['pval_'+args.var1]<0.1).plot(col='time',col_wrap=4,yincrease=False,yscale='log',zorder=0)
dpsis = dtmp[args.var2].where(dtmp['pval_'+args.var2]<0.1)
for a,ax in enumerate(fg.axs.flat):
    dpsis.isel(time=a).plot.contour(levels=vlevels[args.var2],colors='k',ax=ax,x='lat',yincrease=False,yscale='log',zorder=2)
    t = dpsis.isel(time=a).time
    yr = int(t.dt.year.values)
    seas=t.dt.season.values
    if seas == 'DJF': yr += 1
    ax.set_title('year {0} {1}'.format(yr,seas))
    if ax.get_subplotspec().is_first_col():
        ax.set_ylabel('pressure [hPa]')
    else:
        ax.set_ylabel('')
    if ax.get_subplotspec().is_last_row():
        ax.set_xlabel('latitutde')
    else:
        ax.set_xlabel('')
    ax.grid(zorder=3)

ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
outFile = 'figures/{0}_{1}_{2}_QBO_diff_season.pdf'.format(model,args.var1,args.var2)
ac.AddPanelLabels(fg.axs,'upper left',ypos=1.1) 
fg.fig.savefig(outFile,bbox_inches='tight',transparent=True)
print(outFile)
        
