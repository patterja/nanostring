#!/bin/bash
## Usage:
## transfers ftp to Xdrive ../HSR/CLINICAL/Nanostring/automated_data directory

wf_version="1.0.0"


PROJ="${1%/}"
echo "$PROJ"
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

X_DATA_DIR="/Volumes/Histopathology Shared Resource/CLINICAL/Nanostring/automated_data"
X_SS_DIR="/Volumes/Histopathology Shared Resource/CLINICAL/Nanostring/samplesheets"

# X DIRECTORY SYNC
cd "${X_DATA_DIR}"
zfile="${X_DATA_DIR}""/""${PROJ}""_RCC.ZIP"

if [[ -f $zfile ]]; then 
  echo "ZIP already exists."
  exit 1
else
# tabs and spaces seem to mess with heredocs
## Get the data using heredoc
ftp -inv $RAW_HOST > transfer_"${PROJ}".log 2>&1 <<EOF
user $RAW_USER $RAW_PASSWORD
cd RCCData
binary
mget "${PROJ}"_RCC.ZIP
echo "mgetting" "${PROJ}"_RCC.ZIP
ls
bye
EOF

if cat transfer_"${PROJ}".log | grep "The system cannot find the file specified.";  then 
  echo "The system cannot find the file specified. Check project name"
  cat transfer_"${PROJ}".log
else
  ## Decompress data
  unzip "${PROJ}""_RCC.ZIP" -d ${PROJ}
  mv "${PROJ}""_RCC.ZIP" "$HOME/Box Sync/NANOSTRING/RCCData"
  # FIND SAMPLESHEET
  samplesheet=`find "${X_SS_DIR}" -name "${PROJ}*"`
  ss=()
  for i in "$X_SS_DIR"/*;
    do
      case $i in */"${PROJ}"*) 
        echo $i
        ss+=(${i##*/})
        echo ${#ss[@]}
      esac;
  done
  ## check 1 samplesheet and copy to directory
  if [ ${#ss[@]} -eq 0 ]; then
    echo "samplesheet doesn't exist in" 
    echo $X_SS_DIR
    exit 1
  elif [ "${#ss[@]}" -ge 2 ]; then
    echo "multiple samplesheets in Box samplesheet directory, check and rerun"
    exit 1
  else
    cp "$X_SS_DIR""/"$ss "$X_DATA_DIR""/"$PROJ"/"
  fi
fi
fi
#ssh "${REMOTE_SERVER}" "sg clbackup \"rsync -iv -a ${SRC_DIR}/${PROJ}.tar.gz exaclinical:${PROJ_DIR_BASE}\""
