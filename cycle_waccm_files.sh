#!env bash

#$1: directory with netcdf files

source ~/.bashrc
conda activate /g/data/w40/mxj563/software/venv/conda/py3.7
if [ -f hybrid2pressure.template.ncl ]
then
   echo "RE-USING EXISTING hybri2pressure.template.ncl"
else
   cp $repdir/hungatonga/hybrid2pressure.template.ncl .
fi

for file in $1/*.nc
do
    name=${file##*/}
    outFile=${name/.nc/.interp.nc}
    if [ ! -f $outFile ]
    then
       python $repdir/hungatonga/prepare_data.py -i $file -o $outFile
       sh prepare_data.sh
       sh cleanup.sh
    else
	echo "$outFile ALREADY EXISTS. SKIPPING."
    fi
done
