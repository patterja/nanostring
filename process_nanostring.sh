#!/bin/bash

wf_version="1.0.0"


PROJ=${1%/}
echo $PROJ
echo $USER
RAW_HOST="10.125.46.25"
RAW_USER="technician"
RAW_PASSWORD="MAX123"
RAW_DEST="/home/exacloud/clinical/$GROUP/nanostring"
STORAGE_HOST="exanas01.ohsu.edu"
STORAGE_GROUP="clbackup"
STORAGE_DEST_DIR="/home/exacloud/clinicalarchive/KDL/Nanostring"
#REMOTE_HOST="$USER@clinglxy.ohsu.edu"
#REMOTE_DIR="/home/groups/clinical/clinglxy/auto_imports"
REMOTE_HOST="$USER@clinglxydev.ohsu.edu"
REMOTE_DIR="/home/groups/clinical/clinglxydev"

BOX_DATA_DIR="$HOME/Box Sync/NANOSTRING/automated_data"
BOX_SS_DIR="$HOME/Box Sync/NANOSTRING/samplesheets"

# BOX DIRECTORY SYNC
cd "${BOX_DATA_DIR}"
zfile="${BOX_DATA_DIR}""/""${PROJ}""_RCC.ZIP"

if [[ -f $zfile ]]; then 
  echo "ZIP already exists."
  exit 1
else
# tabs and spaces seem to mess with heredocs
## Get the data using heredoc
ftp -inv $RAW_HOST > transfer_${PROJ}.log 2>&1 <<EOF
user $RAW_USER $RAW_PASSWORD
cd RCCData
binary
mget "${PROJ}"_RCC.ZIP
ls
bye
EOF

if cat transfer_${PROJ}.log | grep "The system cannot find the file specified.";  then 
  echo "The system cannot find the file specified. Check project name"
  cat transfer_${PROJ}.log
else
  ## Decompress data
  unzip ${PROJ}"_RCC.ZIP" -d ${PROJ}
  mv ${PROJ}"_RCC.ZIP" "$HOME/Box Sync/NANOSTRING/compressed_data"
  # FIND SAMPLESHEET
  samplesheet=`find "${BOX_SS_DIR}" -name "${PROJ}*"`
  ss=()
  for i in "$BOX_SS_DIR"/*;
    do
      case $i in */${PROJ}*) 
        echo $i
        ss+=(${i##*/})
        echo ${#ss[@]}
      esac;
  done
  ## check 1 samplesheet and copy to directory
  if [ ${#ss[@]} -eq 0 ]; then
    echo "samplesheet doesn't exist in" 
    echo $BOX_SS_DIR
    exit 1
  elif [ ${#ss[@]} -ge 2 ]; then
    echo "multiple samplesheets in Box samplesheet directory, check and rerun"
    exit 1
  else
    cp "$BOX_SS_DIR""/"$ss "$BOX_DATA_DIR""/"$PROJ"/"
  fi
fi

rsync_run()
{
    # Run rsync, put the sequence data in a backed up location.
    timestamp "rsyncing to exaclinical"
    if sg $STORAGE_GROUP -c '/usr/local/bin/rsync -variz "${BOX_DATA_DIR}""/""${PROJ}" "${USER}"@$STORAGE_HOST::$STORAGE_DEST_DIR/'; then
        timestamp "rsync finished without error"
    elif [ $? -eq 5 ]; then
        timestamp "Another transfer is currently in process, waiting 5 minutes before trying again."
        sleep 5m
        rsync_run
    else
        echo "ERROR: rsync failure, check logs!"
#| mailx -s "$(hostname) RSYNC ERROR" "${ADMIN_EMAIL}"
        exit 1
    fi
}
#rsync_run

#ssh -q $USER@$STORAGE_HOST << ENDSSH


#rsync -ivza ${PROJ_DI} "${STORAGE_HOST}::ClinicalArchive/KDL/NextSeqMP/{YEAR}/${PROJ}" --exclude ${PROJ_DIR}/${PROJ}

#rsync -variz ${BOX_DATA_DIR}"/"${PROJ}"_RCC" $REMOTE_HOST:${REMOTE_DIR}"/dataset_import_dir/"${USER}"@ohsu.edu"

fi
