#!/bin/bash
#
#SBATCH -J star_salmon_"$1"
#SBATCH -o logs/"%j".out
#SBATCH -e logs/"%j".out

data_root="$1"
data_folder="$2"
out_folder="$3"
./run_STAR_Salmon.py $data_root $data_folder $out_folder

