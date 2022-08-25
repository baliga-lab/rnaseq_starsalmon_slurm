#!/bin/bash
#
#SBATCH -J star_salmon_"$1"
#SBATCH -o logs/"%j".out
#SBATCH -e logs/"%j".out

genome_dir="$1"
data_root="$2"
data_folder="$3"
out_folder="$4"
./run_STAR_Salmon.py $genome_dir $data_root $data_folder $out_folder

