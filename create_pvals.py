import xarray as xr
from aostools import climate as ac

run_parallel=True

#note: KS test not appropriate here (gives weird values) -> problem of statistics? Also, t-test gives almost exactly same result as WC

pert = xr.open_dataset('aqua_sponge_pert_5yr_ens.nc')
ctrl = xr.open_dataset('aqua_sponge_5yr_ens.nc')

pert_zm = pert.mean('lon')
ctrl_zm = ctrl.mean('lon')

pvals = []
pvals_zm = []
for var in ctrl.data_vars:
    print(var)
    if 'pfull' in ctrl[var].coords:
        pval = ac.StatTest(ctrl_zm[var],pert_zm[var],'WC','member',parallel=run_parallel)
        pval = pval.to_dataset()
        pval['var'] = var
        pvals_zm.append(pval)
    else:
        pval = ac.StatTest(ctrl[var],pert[var],'WC','member',parallel=run_parallel)
        pval = pval.to_dataset()
        pval['var'] = var
        pvals.append(pval)
pval_zmx = xr.concat(pvals_zm,'var')
pval_2dx = xr.concat(pvals,'var')

outFile = 'pval_zm.nc'
pval_zmx.to_netcdf(outFile)
print(outFile)
outFile = 'pval_2d.nc'
pval_2dx.to_netcdf(outFile)
print(outFile)
