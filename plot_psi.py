import xarray as xr
import numpy as np
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',default='waccm',help='Choose model/run to plot.')
parser.add_argument('--qbo',dest='qbo',default=None,help='Only take ensemble members which start in given QBO phase. Either + or - if given.')
parser.add_argument('-y',dest='years',default=None,help='average of this range of years. form startYear,stopYear')
args = parser.parse_args()
model = args.model
qmodel = fc.ModName(model)
sns.set_context('paper')
sns.set_style('whitegrid')
colrs = sns.color_palette()

levels = [l for l in np.linspace(-1e10,1e10,21) if l !=0]

sTS = xr.open_dataset('{0}_psis_season_delta_ens.nc'.format(model),decode_times=False).isel(time=slice(1,None))
sCT = xr.open_dataset('{0}_psis_season_ctrl_ens.nc'.format(model),decode_times=False).isel(time=slice(1,None))
sTS,_,_ = fc.CorrectTime(sTS)
sCT,_,_ = fc.CorrectTime(sCT)
if args.qbo is not None:
    ds = xr.open_dataset('{0}_season_ctrl_ens.nc'.format(model),decode_times=False).isel(time=slice(1,None))
    qbo_pos,qbo_neg = fc.CheckQBO(ds,qmodel)
    if args.qbo == '+':
        filtr = qbo_pos
    elif args.qbo == '-':
        filtr = qbo_neg
    sTS = sTS.isel(member=filtr)
    sCT = sCT.isel(member=filtr)


if args.years is None:
    coldim = 'time'
    ttle = ', year '
else:
    years = [int(y) for y in args.years.split(',')]
    yslce = slice('{0:04d}'.format(years[0]),'{0:04d}'.format(years[1]))
    sTS = sTS.sel(time=yslce).groupby('time.season').mean()
    sCT = sCT.sel(time=yslce).groupby('time.season').mean()
    coldim = 'season'
    ttle = ', years {0}-{1}'.format(*years)
    
ncols = 4

pval = ac.StatTest(sTS.psi,0,'T','member',parallel=True)
psid = sTS.psi.mean('member')
psic = sCT.psi.mean('member')


fp = psic.sel(lev=slice(100,None)).plot.contourf(x='lat',yincrease=False,col=coldim,col_wrap=ncols,levels=20)
for a,ax in enumerate(fp.axs.flat):
    psid.isel({coldim:a}).sel(lev=slice(100,None)).plot.contour(levels=levels,x='lat',yincrease=False,ax=ax,colors='k',linewidths=0.5)
    psid.where(pval<0.1).isel({coldim:a}).sel(lev=slice(100,None)).plot.contour(levels=levels,x='lat',yincrease=False,ax=ax,colors='k',linewidths=2)
    if args.years is None:
        yr = psid.isel({coldim:a}).time.dt.year.values
        ss = psid.isel({coldim:a}).time.dt.season.values
        ax.set_title('{0}{1}{2:2d}'.format(ss,ttle,yr))
    else:
        ss = psid.isel({coldim:a}).season.values
        ax.set_title('{0}'.format(ss)+ttle)
    if not ax.get_subplotspec().is_last_row():
        ax.set_xlabel('')
    if not ax.get_subplotspec().is_first_col():
        ax.set_ylabel('')

outFile = 'figures/{0}_psi_season_year.pdf'.format(model)
if args.years is not None:
    outFile = outFile.replace('_year.pdf','_years{0}-{1}.pdf'.format(*years))
if args.qbo is not None:
    outFile = fc.RenameQBOFile(outFile,args.qbo)
#ac.AddPanelLabels(fp.axs,'upper left',ypos=1.2) 
fp.fig.savefig(outFile,bbox_inches='tight',transparent=True)
print(outFile)
        
