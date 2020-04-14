#!/bin/bash

set -e

#cd "/Users/patterja/Box\ Sync/NANOSTRING/data/${smmart}"
NEWBATCH=$(basename "$1")
echo $NEWBATCH
DATA_DIR="/Volumes/OHSU/CLINICAL/Nanostring"

ss_dir="${DATA_DIR}""/""automated_data/""$NEWBATCH"
output_dir="${DATA_DIR}""/""output/""$NEWBATCH"
#files
FILES=("${DATA_DIR}""/""automated_data/"$NEWBATCH"/*RCC")
abfile="${DATA_DIR}""/""REFERENCE_FILES/ANTIBODY_REFERENCE.csv"
known_pos_file="${DATA_DIR}""/""REFERENCE_FILES/knownpositives.txt"
validation3igg_file="${DATA_DIR}""/""REFERENCE_FILES/validation_samples_iggsub_normalized_20200320.txt"
validation2geo_file="${DATA_DIR}""/""REFERENCE_FILES/validation_samples_geosamp_normalized_20200320.txt"
validationraw_file="${DATA_DIR}""/""REFERENCE_FILES/validation_samples_rawdata_20200320.txt"
md_file="/Users/patterja/Box Sync/NANOSTRING/nanostring_metadata.xlsx"
#find samplesheet
echo $NEWBATCH
echo "${FILES[@]}"
ss=$(find "${ss_dir}" -name "${NEWBATCH}*samplesheet.txt")
echo $ss

#if samplesheet exists process
if test -f "$ss"; then
    mkdir -p "${output_dir}" 
    echo "dir made"
    cd "${output_dir}" 
    echo `pwd`
    /Users/patterja/Workspace/nanostring/nanostring/process_rcc.py --samplesheet "$ss" --abfile "${abfile}" ${FILES[@]}
else
    echo "no samplesheet"
    exit 1
fi

if [ -s rawdata.txt ]; then
    /Users/patterja/Workspace/nanostring/nanostring/norm_nanostring.py --rawdata "rawdata.txt" --abfile "${abfile}"
else
    echo "something went wrong with processing rawdata.txt"
    exit 1
fi

echo "finish processing rawdata, moving on to process_batch_qc.R"

if [ -s 2_GEOMEAN_NORMALIZED.tsv ]; then
    /Users/patterja/Workspace/nanostring/nanostring/process_batch_qc.R -i 2_GEOMEAN_NORMALIZED.tsv --pos_file "${known_pos_file}" --validation_file "${validation2geo_file}" --md_file "${md_file}" --ab_ref_file "${abfile}" 
else
    echo "something went wrong with making 3_IGG_SUBTRACTED.tsv"
    exit 1
fi

if [ -s 2_GEOMEAN_NORMALIZED.tsv ]; then
    echo "running batch correction"
#    /Users/patterja/Workspace/nanostring/nanostring/batch_correction_ruv.R -i "rawdata.txt" --validation_file "${validationraw_file}" --md_file "${md_file}" --ab_ref_file "${abfile}"
    /Users/patterja/Workspace/nanostring/nanostring/batch_correction.R -i "rawdata.txt" --validation_file "${validationraw_file}" --md_file "${md_file}" --ab_ref_file "${abfile}"
else
    echo "something went wrong with making 2_GEOMEAN.tsv"
    exit 1
fi



