
source ~/.bashrc
module use /g/data/hh5/public/modules
module load conda/analysis3


cd $arch/work/HungaTonga/

python plot_vertical_SWV_profile.py


cd $arch/work/HungaTonga/analysis/125Tg/waccm/

python $repdir/hungatonga/plot_q.py --qbo a
python $repdir/hungatonga/plot_QBO_diff_season.py -m waccm --v1 Q --v2 psis
python $repdir/hungatonga/plot_tco.py --years 1,4 --ylim -1 4 --qbo -
python $repdir/hungatonga/plot_tco.py --years 1,4 --ylim -1 4 --qbo +
python $repdir/hungatonga/plot_strat_vapor.py -m waccm --qbo -
python $repdir/hungatonga/plot_strat_vapor.py -m waccm --qbo +
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v U -l 300 -L -4 4 17 -s DJF -Y 1,4 --qbo -
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v U -l 300 -L -4 4 17 -s DJF -Y 1,4 --qbo +
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v U -l 300 -L -4 4 17 -s DJF JJA 
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v V -l 300 -L -4 4 17 -s DJF JJA
python $repdir/hungatonga/plot_plevel_season.py -m waccm -v Z -l 300 -L -80 80 17 -s DJF JJA
python $repdir/hungatonga/plot_maps.py -m waccm --center 155 -v TS P
python $repdir/hungatonga/plot_maps.py -m waccm --center 155 -v DLS CLDTOT LWCF SWCF
python $repdir/hungatonga/plot_ipv_waf.py -m waccm -s DJF JJA -y 3,7
python $repdir/hungatonga/plot_lines_regions.py -m waccm -v TS

cd ../mima
python $repdir/hungatonga/plot_maps.py -m bench_SH_fix_125Tg --center 155 -v TS
python $repdir/hungatonga/plot_maps.py -m bench_SH_fixq_o3init_125Tg --center 155 -v TS
python $repdir/hungatonga/plot_maps.py -m bench_SH_fixo3_qinit_125Tg --center 155 -v TS
