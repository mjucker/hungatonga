
source ~/.bashrc
module use /g/data/hh5/public/modules
module load conda/analysis3


cd $arch/work/HungaTonga/
# Figure 1 
python $repidr/hungatonga/plot_vertical_SWV_profile.py

cd $arch/work/HungaTonga/analysis/125Tg/waccm/

# Figure 2
python $repdir/hungatonga/plot_q.py --qbo a
# Figure 3
python $repdir/hungatonga/plot_QBO_diff_season.py -m waccm --v1 Q --v2 psis
# Figure 5ab
python $repdir/hungatonga/plot_tco.py --years 1,4 --ylim -1 4 --qbo -
python $repdir/hungatonga/plot_tco.py --years 1,4 --ylim -1 4 --qbo +
# Figure 5cd
python $repdir/hungatonga/plot_strat_vapor.py -m waccm --qbo -
python $repdir/hungatonga/plot_strat_vapor.py -m waccm --qbo +
# additional figure to compare MLS and WACCM 2023-10
python $repdir/hungatonga/plot_strat_vapor.py -m waccm --qbo - -L 0001-12
python $repdir/hungatonga/plot_strat_vapor.py -m waccm --qbo - -L 0001-10


python $repdir/hungatonga/plot_plevel_season.py -m waccm -v U -l 300 -L -4 4 17 -s DJF -Y 1,4 --qbo -
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v U -l 300 -L -4 4 17 -s DJF -Y 1,4 --qbo +
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
# Figure B2
python $repdir/hungatonga/plot_lines_regions.py -m waccm -v TS
python $repdir/hungatonga/plot_zm_season.py -m waccm -s DJF JJA -y 3,7 -v T -L -0.5 0.5 26
python $repdir/hungatonga/plot_zm_season.py -m waccm -s DJF JJA -y 3,7 -v U
# extra figure: global mean temperature anomalies
python $repdir/hungatonga/compute_global_mean_temperature.py -m waccm
# extra figure: bootstrapped means
python $repdir/hungatonga/bootstrap_members.py -m waccm -v TS
python $repdir/hungatonga/bootstrap_members.py -m waccm -v P
# extra figure: TS every year
python $repdir/hungatonga/plot_maps_years.py -m waccm --center 155 -v TS -y 0 9 --seasons DJF
python $repdir/hungatonga/plot_maps_years.py -m waccm --center 155 -v TS -y 0 9 --seasons JJA
# extra figure: PV and WAF every year
#python $repdir/hungatonga/plot_ipv_waf.py -m waccm -s DJF -y 0-9 -l 0 --sig
#python $repdir/hungatonga/plot_ipv_waf.py -m waccm -s JJA -y 0-9 -l 0 --sig



## MiMA plots
cd ../mima
# Figure 11
python $repdir/hungatonga/plot_maps.py -m bench_SH_fix_125Tg --center 155 -v TS P -b
# supplementary figures
python $repdir/hungatonga/plot_plevel_season.py -m bench_SH_fix_125Tg -v CFSD -l 200 -y 3,7 -s DJF JJA -C PRGn_r -L -0.6 0.6 10
python $repdir/hungatonga/plot_maps.py -m bench_SH_fixq_o3init_125Tg --center 155 -v TS P -b
python $repdir/hungatonga/plot_maps.py -m bench_SH_fixo3_qinit_125Tg --center 155 -v TS P -b
python $repdir/hungatonga/plot_plevel_season.py -m bench_SH_fix_125Tg -v V -l 300 -L -3 3 16 -y 3,7
# global mean temperature anomalies
python $repdir/hungatonga/compute_global_mean_temperature.py -m bench_SH_fix_125Tg


