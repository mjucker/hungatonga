#!env bash

# this file will assume the order of restart_00??.tar corresponds
#  one-to-one with the order of the member files
# Thus, to make this work, first copy all necessary restart_00??.tar files
#  into the working directory

basedir=$1

suffix=$2


member_dir=${basedir}/members/pert

if [ -d RESTART ]
then
    rm -rf RESTART/*
else
    mkdir RESTART
fi

tars=(restart_00??.tar)
q_membs=("${member_dir}"/*_q_*)
o_membs=("${member_dir}"/*_o3_*)
n=0
for file in "${tars[@]}"
do
    cp "${q_membs[n]}" RESTART/q${suffix}.nc
    cp "${o_membs[n]}" RESTART/o3${suffix}.nc
    tar rf $file RESTART/*
    rm RESTART/*
    echo $file
    let n=$n+1
done
