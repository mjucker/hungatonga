import xarray as xr
from hungatonga import functions as fc
from glob import glob

fles = glob('*_season_delta_ens.nc')
fles.sort()


for fle in fles:
    if 'cumflx' in fle:
        continue
    ds = xr.open_dataset(fle,decode_times=False)
    ds,_,_ = fc.CorrectTime(ds)
    model = fle.replace('_season_delta_ens.nc','')
    cdv = []
    for var in ['fsd','fld']:
        cds = []
        for p,pres in enumerate(ds.pfull.values):
            if p == 0:
                cfsd = ds[var].sel(pfull=pres)
            else:
                cfsd = -ds[var].sel(pfull=slice(pres,None)).integrate('pfull')/pres
                cfsd['pfull'] = pres
            cfsd.name = 'c'+cfsd.name
            cds.append(cfsd)
        cdv.append(xr.concat(cds,'pfull'))
    dt = xr.merge(cdv)
    outFile = model+'_cumflx_season_delta_ens.nc'
    dt.to_netcdf(outFile)
    print(outFile)
