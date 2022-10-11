#!env bash

#$1: directory with netcdf files

source ~/.bashrc
conda activate /g/data/w40/mxj563/software/venv/conda/py3.7
cp $repdir/hungatonga/hybrid2pressure.template.ncl .

for file in $1/*.nc
do
    name=${file##*/}
    outFile=${name/.nc/.interp.nc}
    python $repdir/hungatonga/prepare_data.py -i $file -o $outFile
    sh prepare_data.sh
    sh cleanup.sh
done
