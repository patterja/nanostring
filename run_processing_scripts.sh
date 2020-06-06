#!/bin/bash
#
# Pipeline to run everything together. 
## Usage as loop
## for batch in /Volumes/OHSU/CLINICAL/Nanostring/automated_data/*; do echo $batch; run_processing_scripts.sh $batch; done
##
## Usage:
## run_processing_scripts.sh 20190314_208421591020
##
set -e

valid_version="20200512"
ihc_version="20200507"
knownpos_version="1.0"
antibody_version="1.0"

#cd "/Users/patterja/Box\ Sync/NANOSTRING/data/${smmart}"
NEWBATCH=$(basename "$1")
echo $NEWBATCH
DATA_DIR="/Volumes/OHSU/CLINICAL/Nanostring"
REPO_DIR="/Users/patterja/Workspace/nanostring/nanostring"

ss_dir="${DATA_DIR}""/""automated_data/""$NEWBATCH"
output_dir="${DATA_DIR}""/""output/""$NEWBATCH"
#files
FILES=("${DATA_DIR}""/""automated_data/"$NEWBATCH"/*RCC")
abfile="${DATA_DIR}""/""REFERENCE_FILES/ANTIBODY_REFERENCE_v"${antibody_version}".csv"
known_pos_file="${DATA_DIR}""/""REFERENCE_FILES/knownpositives_v"${knownpos_version}".txt"
#validation3igg_fle="${DATA_DIR}""/""REFERENCE_FILES/validation_samples_iggsub_normalized_20200320.txt"
#validation2geo_file="${DATA_DIR}""/""REFERENCE_FILES/validation_samples_geosamp_normalized_"${valid_version}".txt"
validationraw_file="${DATA_DIR}""/""REFERENCE_FILES/validation_samples_rawdata_"${valid_version}".txt"
ihc_file="${DATA_DIR}""/""REFERENCE_FILES/ihc_status_"${ihc_version}".txt"
md_file="/Users/patterja/Box Sync/NANOSTRING/nanostring_metadata.xlsx"
#find samplesheet
echo $NEWBATCH
echo "${FILES[@]}"
ss=$(find "${ss_dir}" -name "${NEWBATCH}*samplesheet.txt")
echo $ss

#Step 1: if samplesheet exists process raw RCC files
if test -f "$ss"; then
    mkdir -p "${output_dir}" 
    echo "dir made"
    cd "${output_dir}" 
    echo `pwd`
    python3 "${REPO_DIR}""/"process_rcc.py --samplesheet "$ss" --abfile "${abfile}" ${FILES[@]}
else
    echo "no samplesheet"
    exit 1
fi

#Step 2: if rawdata.txt exists run per batch scaling 
if [ -s rawdata.txt ]; then
    python3 "${REPO_DIR}""/"norm_nanostring.py --rawdata "rawdata.txt" --abfile "${abfile}"
else
    echo "something went wrong with processing rawdata.txt"
    exit 1
fi

echo "finish processing rawdata, moving on to process_batch_qc.R"

#Step 3: if rawdata.txt exists process QC stats for batch pass/fail
if [ -s rawdata.txt ]; then
    "${REPO_DIR}""/"process_batch_qc.R -i "rawdata.txt" --pos_file "${known_pos_file}" --validation_file "${validationraw_file}" --md_file "${md_file}" --ab_ref_file "${abfile}" 
else
    echo "something went wrong with making rawdata.txt"
    exit 1
fi

#Step 4: if rawdata.txt exists process sample batch correction with cohort
if [ -s rawdata.txt ]; then
    echo "running batch correction"
#    /Users/patterja/Workspace/nanostring/nanostring/batch_correction_ruv.R -i "rawdata.txt" --validation_file "${validationraw_file}" --md_file "${md_file}" --ab_ref_file "${abfile}"
    "${REPO_DIR}""/"batch_correction_ruv.R --input "rawdata.txt" --validation_file "${validationraw_file}" --md_file "${md_file}" --ihc_file "${ihc_file}" --ab_ref_file "${abfile}"
else
    echo "something went wrong with making 2_GEOMEAN.tsv"
    exit 1
fi



