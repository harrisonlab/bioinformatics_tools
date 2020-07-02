#!/usr/bin/env bash
#SBATCH -J cassis
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4



Assembly=$1
CassisTSV=$2
GeneID=$3
OutDir=$4

cassis \
  --annotation $CassisTSV \
  --genome $Assembly \
  --anchor $GeneID \
  --dir $OutDir/$GeneID \
  --mismatches 0 \
  -v \
  --prediction \
  --num-cpus 4 \
  | tee 2>&1 $OutDir/$GeneID/${GeneID}_log.txt
