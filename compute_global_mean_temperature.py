import xarray as xr
from hungatonga import functions as fc
from aostools import climate as ac
import argparse
from matplotlib import pyplot as plt
import seaborn as sns
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',help='Choose model/run to plot.')
parser.add_argument('-y',dest='years',default='3,7',help='average of this range of years. form startYear,stopYear')
args = parser.parse_args()
model = args.model
qmodel = fc.ModName(model)
sns.set_context('notebook')
sns.set_style('whitegrid')
colrs = sns.color_palette()

du = xr.open_dataset(args.model+'_season_delta_ens.nc',decode_times=False)
ds,_,_ = fc.CorrectTime(du)

years = [int(y) for y in args.years.split(',')]
yslce = slice('{0:04d}'.format(years[0]),'{0:04d}'.format(years[1]))

if 'TS' in ds.data_vars:
    var = 'TS'
else:
    var = fc.variables[qmodel]['TS']

tg = ac.GlobalAvgXr(ds[var].mean('lon'))
th = ac.GlobalAvgXr(du[var].mean('lon'))

ta = tg.sel(time=yslce).mean('time')
std = ta.std('member')
mean= ta.mean('member')

print('Global mean Ts anomaly for years {0}-{1} is {2:5.3f}+-{3:5.3f}K'.format(years[0],years[1],mean,std))


th = th.rolling(time=3,center=True).mean()
if 'waccm' in qmodel:
    daysperyear = 365
else:
    daysperyear = 360
th = th.assign_coords({'time':th.time/daysperyear})
th.time.attrs['units'] = 'years since eruption'
pval = ac.StatTest(th,0,'T','member')
tstd = th.std('member')
tmean= th.mean('member')
fig,ax = plt.subplots()
ax.axhline(0,color='k')
p = tmean.plot(ax=ax,lw=1,ls='--')[0]
tmean.where(pval<0.1).plot(ax=ax,color=colrs[0])
ax.fill_between(p.get_xdata(),tmean-tstd,tmean+tstd,color=colrs[0],alpha=0.3)
ax.set_xlim(0,10)
ax.set_ylabel('global mean temperature anomalies [K]')
ax.set_title('{0} global mean Ts anomalies, 3-month rolling mean [K]'.format(qmodel.upper()))
outFile = 'figures/{0}_global_mean_ts.pdf'.format(args.model)
fig.savefig(outFile,transparent=True,bbox_inches='tight')
print(outFile)


