#!/bin/bash

# Submit the R directories in the specified file to the cluster

if [ $# -lt 3 ]
then
  echo "Usage: submit_list.sh <data_root> <rnum_list> <outfolder>"
  exit 0
fi

DATA_ROOT=$1
OUT_FOLDER=$3

for d in `cat $2` ; do
  data_folder="$d"
  echo "submitting $data_folder"
  sbatch ./run_star_salmon_job.sh $DATA_ROOT $data_folder $OUT_FOLDER
done
