#env bash
module load nco

rm -rf RESTART
for memb in {$1..$2}
do
    rm out.00${memb}.txt
    tar xf restart_00${memb}.tar
    ncbo -y add RESTART/spectral_dynamics.res.nc delta_init.nc -O RESTART/spectral_dynamics.res.nc
    tar cf restart_00${memb}.tar
    echo $memb
done

