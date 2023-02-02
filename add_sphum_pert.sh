#env bash
set -e
module load nco

rm -rf RESTART
for ((memb=$1; memb<=$2; memb++))
do
    #outfile=out.00${memb}.txt
    #if [[ -f $outfile ]]
    #then
#	rm $outfile
    #fi
    tar xf restart_00${memb}.tar
    ncbo -y add RESTART/spectral_dynamics.res.nc delta_init.nc -O RESTART/spectral_dynamics.res.nc
    tar cf restart_00${memb}.tar RESTART
    echo $memb
done

