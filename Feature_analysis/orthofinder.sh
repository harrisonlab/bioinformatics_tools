#!/usr/bin/env bash
#SBATCH -J orthofinder
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=24

# Find orthogroups and orthologs


IN_DIR=$1
prefix=$2

cd $IN_DIR

orthofinder -f ./ -t 24 -a 6 -n $prefix 