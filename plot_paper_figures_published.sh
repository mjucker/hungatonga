repdir=<path/to/hungatonga_repository>

data_dir=<path/to/data/dir>

cd $data_dir
# Figure 1 
python $repdir/hungatonga/plot_vertical_SWV_profile.py

# Figures 2,4
python $repdir/hungatonga/plot_q.py --qbo a
# Figure 3
python $repdir/hungatonga/plot_QBO_diff_season.py -m waccm --v1 Q --v2 psis
# Figure 5ab
python $repdir/hungatonga/plot_tco.py --years 1,4 --ylim -1 4 --qbo -
python $repdir/hungatonga/plot_tco.py --years 1,4 --ylim -1 4 --qbo +
# Figure 5cd
python $repdir/hungatonga/plot_strat_vapor.py -m waccm --qbo -
python $repdir/hungatonga/plot_strat_vapor.py -m waccm --qbo +
# Supplementary Figure S5
python $repdir/hungatonga/plot_strat_vapor.py -m waccm --qbo - -L 0001-12
python $repdir/hungatonga/plot_strat_vapor.py -m waccm --qbo - -L 0001-10
# Figure 6abcd
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v U -l 300 -L -4 4 17 -s DJF -Y 1,4 --qbo -
# Figure 6efgh
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v U -l 300 -L -4 4 17 -s DJF -Y 1,4 --qbo +
# Figure S9
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v U -l 300 -L -4 4 17 -s DJF JJA -y 3,7
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v V -l 300 -L -4 4 17 -s DJF JJA -y 3,7
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v Z -l 300 -L -80 80 17 -s DJF JJA -y 3,7
# Figures 7,8
python $repdir/hungatonga/plot_maps.py -m waccm --center 155 -v TS P
# Figure 9a-d
python $repdir/hungatonga/plot_maps.py -m waccm --center 155 -v DLSC TCWV
# Figure 9ef
python $repdir/hungatonga/plot_ipv_waf.py -m waccm -s DJF JJA -y 3,7 -l 4
# Figure 10
python $repdir/hungatonga/plot_maps.py -m waccm --center 155 -v CLDTOT LWCF SWCF
# Figure S6
python $repdir/hungatonga/plot_lines_regions.py -m waccm -v TS
# Figure S7
python $repdir/hungatonga/plot_maps_years.py -m waccm --center 155 -v TS -y 0 9 --seasons DJF --sig
# Figure S8
python $repdir/hungatonga/plot_maps_years.py -m waccm --center 155 -v TS -y 0 9 --seasons JJA --sig

## MiMA plots
# Figure 11
python $repdir/hungatonga/plot_maps.py -m bench_SH_fix_125Tg --center 155 -v TS P -b
# Figure S13
python $repdir/hungatonga/plot_plevel_season.py -m bench_SH_fix_125Tg -v CFSD -l 200 -y 3,7 -s DJF JJA -C PRGn_r -L -0.6 0.6 10
# Figure S11
python $repdir/hungatonga/plot_maps.py -m bench_SH_fixq_o3init_125Tg --center 155 -v TS P -b
# Figure S12
python $repdir/hungatonga/plot_maps.py -m bench_SH_fixo3_qinit_125Tg --center 155 -v TS P -b
# Figure S10
python $repdir/hungatonga/plot_plevel_season.py -m bench_SH_fix_125Tg -v V -l 300 -L -3 3 16 -y 3,7


