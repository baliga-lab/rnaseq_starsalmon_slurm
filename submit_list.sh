#!/bin/bash

# Submit the R directories in the specified file to the cluster

if [ $# -lt 4 ]
then
  echo "Usage: submit_list.sh <genome_dir> <data_root> <rnum_list> <outfolder>"
  exit 0
fi

GENOME_DIR=$1
DATA_ROOT=$2
OUT_FOLDER=$4

for d in `cat $3` ; do
  data_folder="$d"
  echo "submitting $data_folder"
  sbatch ./run_star_salmon_job.sh $GENOME_DIR $DATA_ROOT $data_folder $OUT_FOLDER
done
