import xarray as xr
import numpy as np
from hungatonga import functions as fc
from aostools import climate as ac
from hungatonga.plot_maps_months import vmaxs,labls,cmaps
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
parser.add_argument('-v',dest='variable',default='TS',help='Choose variable to plot.')
parser.add_argument('-M',dest='nmembs',default=15,type=int,help='Set number of members per bootstrap iteration.')
parser.add_argument('-n',dest='nboots',default=1000,type=int,help='Set number of bootstrap members.')


args = parser.parse_args()
model = fc.ModName(args.model)
var = fc.variables[model][args.variable]
nboots = args.nboots

scales = {}
if 'waccm' in args.model.lower():
    scales['P'] = 1e3*86400
    scales['SLP'] = 0.01
else:
    scales['P'] = 86400
    scales['SLP'] = 1

np.random.seed(42)


ds = xr.open_dataset(model+'_season_delta_ens.nc',decode_times=False)
ds,_,_ = fc.CorrectTime(ds)

dTS = ds[var].sel(time=slice('0003','0007')).groupby('time.season').mean()
nmembs = len(dTS.member)
if args.variable in scales.keys():
    dTS =dTS*scales[args.variable]

boots = []
for n in range(nboots):
    sel = np.random.choice(np.arange(nmembs),args.nmembs)
    boot = dTS.isel(member=sel).mean('member')
    boot['n'] = n
    boots.append(boot)
bTS = xr.concat(boots,'n')
pval = ac.StatTest(bTS,0,'T','n',parallel=True)
bTS = bTS.mean('n')

fig,axs,transf = ac.Projection('PlateCarree',ncols=2,nrows=2,kw_args={'central_longitude':155})
if args.variable in cmaps.keys():
    transf['cmap'] = cmaps[args.variable]
fig.set_figwidth(12)
for a,ax in enumerate(axs.flat):
    tmp = bTS.isel(season=a)
    seas = tmp.season.values
    cf=tmp.plot(ax=ax,vmax=vmaxs[args.variable],add_colorbar=False,**transf)
    ax.set_title(seas)
    ax.gridlines()
    ax.coastlines()
ac.AddColorbar(fig,axs,cf,shrink=0.95,cbar_args={'label':labls[args.variable]})
outFile = 'figures/{0}_{1}_bootstrap_{2}_{3}.pdf'.format(args.model,args.variable,args.nmembs,args.nboots)
fig.savefig(outFile,bbox_inches='tight',transparent=True)
print(outFile)
    
