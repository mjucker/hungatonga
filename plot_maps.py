#!env python
import xarray as xr
import numpy as np
import seaborn as sns
from aostools import climate as ac
from hungatonga import functions as fc
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator,FormatStrFormatter
from shapely import geometry
from cartopy import crs as ccrs
import calendar
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',dest='model',help='Choose model/run to plot.')
parser.add_argument('--qbo',dest='qbo',default=None,help='Only take ensemble members which start in given QBO phase. Either + or - if given.')
parser.add_argument('--qmodel',dest='qbo_model',default=None,help='Use this model to assess QBO phase.')
parser.add_argument('--seasons',dest='season_overwrite',default=None,nargs='+',help='Overwrite default seasons to plot.')
parser.add_argument('--center',dest='center_overwrite',default=None,type=float,help='Overwrite default longitude center of plot.')
parser.add_argument('-v',dest='vars',default=None,nargs='+',help='Plot these variables.')
parser.add_argument('-b',dest='boxes',action='store_false',help='Do NOT plot regional boxes in the maps, and do NOT plot bar plots.')
args = parser.parse_args()
model = args.model
qmodel = fc.ModName(model)
sns.set_context('notebook')
sns.set_style('whitegrid')
colrs = sns.color_palette()

pthresh = 0.1

if args.vars is None:
    fields = ['TS','P','SLP','OLR','DLS','CLDTOT']
else:
    fields = args.vars

trans_fields = {'TS':'TREFHT'}

scales = {'TS': 1, 'OLR': 1, 'DLS': 1, 'CLDTOT': 1, 'LWCF': 1, 'SWCF': 1, 'ICEFRAC': 1,'DSS': 1}
for var in fields:
    if 'CLD' in var:
        scales[var] = scales['CLDTOT']

if 'waccm' in model:
    scales['P'] = 1e3*86400
    scales['SLP'] = 0.01
else:
    scales['P'] = 86400
    scales['SLP'] = 1

default_seasons = ['DJF','MAM','JJA','SON']
seasons = {
    'TS'    :default_seasons,
    'OLR'   :['DJF','JJA'],
    'DLS'   :['DJF','JJA'],
    'DSS'   :['DJF','JJA'],
    'P'     :default_seasons,
    'CLDTOT':['DJF','JJA'],
    'LWCF'  :['DJF','JJA'],
    'SWCF'  :['DJF','JJA'],
    'SLP'   :default_seasons,
#    'SLP'   :[1,7],
    'ICEFRAC':default_seasons,
}
if args.season_overwrite is not None:
    aseas = args.season_overwrite
    for key in seasons.keys():
        seasons[key] = aseas
for var in fields:
    if 'CLD' in var:
        seasons[var] = seasons['CLDTOT']

nseas = 0
for field in fields:
    nseas = max(nseas,len(seasons[field]))

cent_def = 0 #155
proj = 'PlateCarree' #'Robinson'

loncents = {
    'P' : cent_def,
    'OLR':155, #cent_def,
    'DLS':155, #cent_def,
    'DSS':155, #cent_def,
    'LWCF':155, #cent_def,
    'SWCF':155, #cent_def,
    'TS': cent_def,
    'CLDTOT':155,#cent_def,
    'SLP':cent_def,
    'ICEFRAC':cent_def,
}
if args.center_overwrite is not None:
    for key in loncents.keys():
        loncents[key] = args.center_overwrite
for var in fields:
    if 'CLD' in var:
        loncents[var] = loncents['CLDTOT']

def_yrs = [3,7]
avg_years = {
    'P':   def_yrs,
    'OLR': def_yrs,
    'DLS': def_yrs,
    'DSS': def_yrs,
    'TS':  def_yrs,
    'CLDTOT':  def_yrs,
    'LWCF':  def_yrs,
    'SWCF':  def_yrs,
    'SLP': def_yrs,
    'ICEFRAC': def_yrs,
}
for var in fields:
    if 'CLD' in var:
        avg_years[var] = avg_years['CLDTOT']

vmaxs = {
    'TS': 1.5,
    'P' : 1.0,
    'SLP': 6,
    'CLDTOT': 0.05,
    'OLR': 5.0,
    'DLS': 5.0,
    'DSS': 5.0,
    'LWCF': 5.0,
    'SWCF': 5.0,
    'ICEFRAC': None,
    }
for var in fields:
    if 'CLD' in var:
        vmaxs[var] = vmaxs['CLDTOT']

cmaps = {
    'TS': 'RdBu_r',
    'P' : 'PuOr',
    'OLR':'PiYG',
    'DLS':'PRGn_r',
    'DSS':'PRGn_r',
    'LWCF':'PRGn_r',#'cool',
    'SWCF':'PRGn_r',#'cool',
    #'CLDTOT':'ocean',
    'CLDTOT':'cividis',
    'SLP':'BrBG_r',
    'ICEFRAC':'cubehelix',
    }
for var in fields:
    if 'CLD' in var:
        cmaps[var] = cmaps['CLDTOT']

labls = {
    'TS': 'T2m [K]',
    'P' : 'Q [mm/day]',
    'OLR':'OLR [W/m2]',
    'DLS':'DLS [W/m2]',
    'DSS':'DSS [W/m2]',
    'LWCF':'LWCF [W/m2]',
    'SWCF':'SWCF [W/m2]',
    'CLDTOT':'CLD []',
    'SLP':'SLP [hPa]',
    'ICEFRAC':'[]',
    }
for var in fields:
    if 'CLD' in var:
        labls[var] = labls['CLDTOT']

areas = fc.areas

#pert = xr.open_dataset('{0}_pert_ens.nc'.format(model),decode_times=False)
#ctrl = xr.open_dataset('{0}_ctrl_ens.nc'.format(model),decode_times=False)
if args.qbo is not None:
    pertq = xr.open_dataset(args.qbo_model+'_pert_ens.nc',decode_times=False)
    if args.qbo_model is None:
        qbo_pos,qbo_neg = fc.CheckQBO(pertq,model)
    else:
        qbo_pos,qbo_neg = fc.CheckQBO(pertq,args.qbo_model)
        qbo_pos = qbo_pos.assign_coords({'member':pert.member})
        qbo_neg = qbo_neg.assign_coords({'member':pert.member})
    if args.qbo == '+':
        pert = pert.isel(member=qbo_pos)
        ctrl = ctrl.isel(member=qbo_pos)
    elif args.qbo == '-':
        pert = pert.isel(member=qbo_neg)
        ctrl = ctrl.isel(member=qbo_neg)


#dTS = pert - ctrl
dTS = xr.open_dataset(args.model+'_season_delta_ens.nc',decode_times=False)
dTS,_,_ = fc.CorrectTime(dTS)
end_year = '{0:04d}'.format(dTS.time[0].dt.year+9)
dTS = dTS.sel(time=slice(None,end_year))

nvars = len(fields)
#ncols = nseas
ncols = 2
#nrows = nvars
nrows = int(nseas/ncols*nvars)
fig = plt.figure(figsize=[6*ncols,3*nrows])
transf = {'transform':ccrs.PlateCarree()}

axs = []
p = 0
for f,field in enumerate(fields):
    var = fc.variables[qmodel][field]
    if var in trans_fields.keys():
        var = trans_fields[var]
    da = dTS[var]*scales[field]
    # first, seasonal maps
    proj_func = getattr(ccrs,proj)
    faxs = []
    for s,season in enumerate(seasons[field]):
        p += 1
        ax = fig.add_subplot(nrows,ncols,p,projection=proj_func(central_longitude=loncents[field]))
        axs.append(ax)
        faxs.append(ax)
        if isinstance(season,int):
            filtr = da['time.month'] == season
            seas = calendar.month_abbr[season]
        else:
            filtr = da['time.season'] == season
            seas = season
        dm = da.isel(time=filtr)
        dm = dm.sel(time=slice('{0:04d}'.format(avg_years[field][0]),'{0:04d}'.format(avg_years[field][1]))).mean('time')
        pval = ac.StatTest(dm,0,'T','member',parallel=True)
        dmm = dm.mean('member')
        cf = dmm.where(pval<pthresh).plot(ax=ax,vmax=vmaxs[field],cmap=cmaps[field],add_colorbar=False,**transf)
        ax.gridlines()
        ax.coastlines()
        if args.boxes and field in areas.keys():
            geoms = []
            r=0
            for region,sels in areas[field][season].items():
                if sels['lon'].start > 180 and sels['lon'].stop-sels['lon'].start < 180:
                    lonstart = sels['lon'].start - 360
                else:
                    lonstart = sels['lon'].start
                if sels['lon'].stop > 180 and sels['lon'].stop-sels['lon'].start < 180:
                    lonstop = sels['lon'].stop - 360
                else:
                    lonstop = sels['lon'].stop
                geom = geometry.box(minx=lonstart,maxx=lonstop,
                                    miny=sels['lat'].start,maxy=sels['lat'].stop)
                #geoms.append(geom)
                ax.add_geometries([geom],edgecolor=colrs[r],facecolor='none',crs=ccrs.PlateCarree(),zorder=5)
                r += 1
        ax.set_title(field+' anomalies, {0}, years {1}-{2}'.format(seas,avg_years[field][0],avg_years[field][1]))
    ac.AddColorbar(fig,faxs,cf,shrink=0.85,cbar_args={'label':labls[field]})
#if ncols == nseas:
ac.AddPanelLabels(axs,'upper left',ypos=1.15) 
#else:
#    ac.AddPanelLabels(axs,xpos=-0.1,ypos=0.1) 
outFile = 'figures/{0}_maps_boxed.pdf'.format(model)
if args.vars is not None:
    outFile = outFile.replace('.pdf','_'+'_'.join(fields)+'.pdf')
outFile = outFile.replace('.pdf','_{0:d}-{1:d}.pdf'.format(*avg_years[field]))
if args.qbo is not None:
    outFile = fc.RenameQBOFile(outFile,args.qbo)
fig.savefig(outFile,bbox_inches='tight',transparent=True)
print(outFile)

# now add bar charts for the regions
if args.boxes:
    for field in fields:
        if field not in areas.keys():
            continue
        var = fc.variables[qmodel][field]
        if var in trans_fields.keys():
            var = trans_fields[var]
        da = dTS[var]*scales[field]
        for season in areas[field].keys():
            dt = []
            nbars = 0
            filtr = da['time.season'] == season
            das = da.isel(time=filtr)
            for region,sels in areas[field][season].items():
                tmp = ac.GlobalAvgXr(das.sel(sels)).mean('lon')
                labl = ' '.join([region,season])
                tmp = tmp.groupby('time.year').mean()
                tmp['region'] = labl
                dt.append(tmp)
                nbars += 1
            if nbars == 0:
                continue
            fig,ax = plt.subplots()
            dx = xr.concat(dt,'region')
            dx.name = field
            df = fc.MakeDataFrame(dx,'region')
            fig.set_figwidth(4*(1+0.5*(nbars-1)))
            try: #new versions of sns
                sx = sns.barplot(data=df,x='year',y=field,hue='region',errorbar=('ci',90),ax=ax)
            except:
                sx = sns.barplot(data=df,x='year',y=field,hue='region',ci=90,ax=ax)
            sns.move_legend(sx,[1.01,0.05])
            #ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            #ax.grid(True,axis='x',which='minor')
            nyrs = len(dx.year)
            xticks = np.linspace(*ax.get_xlim(),nyrs+1)
            ax.set_xticks(xticks)
            xlabels = ['{0:g}'.format(np.abs(np.ceil(l))) for l in xticks]
            ax.set_xticklabels(xlabels)
            ax.set_xlim(xticks[0],xticks[-1])
            ax.grid(True,axis='x')
            sns.despine(bottom=True,left=True)
            ax.set_title('{0} anomalies, {1}'.format(field,season))
            ax.set_ylabel(labls[field])
            ax.set_ylim(-vmaxs[field]*1.5,vmaxs[field]*1.5)
            ax.set_xlabel('time [years since eruption]')
            ylims = ax.get_ylim()
            ax.plot([avg_years[field][0]-0.5,avg_years[field][1]+0.5],[ylims[0],ylims[0]],color='r',lw=2)
            outFile = 'figures/{0}_bars_{1}_{2}.pdf'.format(model,field,season)
            if args.qbo is not None:
                outFile = fc.RenameQBOFile(outFile,args.qbo)
            fig.savefig(outFile,bbox_inches='tight',transparent=True)
            print(outFile)


