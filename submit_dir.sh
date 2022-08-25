#!/bin/bash

# Submit the specified data root subdirectories to the cluster

if [ $# -lt 3 ]
then
  echo "usage submit_dir.sh <genome_dir> <data_root> <outfolder>"
  exit 0
fi

GENOME_DIR=$1
DATA_ROOT=$2
OUT_FOLDER=$3
#"/proj/omics4tb2/Global_Search/Pilot_Pass/X204SC21081158-Z01-F003/raw_data"

for d in "$DATA_ROOT"/*/ ; do
  if [ -d "$d" ]
  then
    data_folder=`basename $d`
    echo "submitting $data_folder"
    sbatch ./run_star_salmon_job.sh $GENOME_DIR $DATA_ROOT $data_folder $OUT_FOLDER
  else
    echo "Skipping $data_folder (not a directory)"
  fi
done
