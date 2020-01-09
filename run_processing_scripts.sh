#!/bin/bash

#cd "/Users/patterja/Box\ Sync/NANOSTRING/data/${smmart}"
NEWBATCH=`basename $1`
DATA_DIR="/Volumes/Histopathology Shared Resource/CLINICAL/Nanostring"

abfile="${DATA_DIR}""/""REFERENCE_FILES/ANTIBODY_REFERENCE.csv"
ss_dir="${DATA_DIR}""/""automated_data/""$NEWBATCH"
output_dir="${DATA_DIR}""/""output_subtracted_newfig/""$NEWBATCH"
FILES=("${DATA_DIR}""/""automated_data/""$NEWBATCH"/*RCC)


#find samplesheet
echo $NEWBATCH
echo "${FILES[@]}"
ss=`find "$ss_dir" -name "${NEWBATCH}*samplesheet.txt"`
echo $ss

#if samplesheet exists process
if test -f "$ss"; then
    mkdir -p "${output_dir}" 
    echo "dir made"
    cd "${output_dir}" 
    /Users/patterja/Workspace/nanostring/nanostring/proc_norm.py --samplesheet "$ss" --abfile "${abfile}" "${FILES[@]}"
else
    echo "no samplesheet"
fi

if [ -s 3_IGG_SUBTRACTED.tsv ]; then
    /Users/patterja/Workspace/nanostring/nanostring/ruv_mbc.R -i 3_IGG_SUBTRACTED.tsv
else
    echo "something went wrong with making 3_IGG_SUBTRACTED.tsv"
fi
