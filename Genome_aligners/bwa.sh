#!/usr/bin/env bash
#SBATCH -J bwa
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

Prefix=$1
Assembly=$2
F_reads=$3
R_reads=$4
OutDir=$5

CWD=$PWD

WorkDir="$TMPDIR"
mkdir -p $WorkDir

cp -r $Assembly $F_reads $R_reads $WorkDir
Fr=$(basename "$F_reads")
Rr=$(basename "$R_reads")
As=$(basename "$Assembly")

cd $WorkDir

# Align the data
bwa index $As

bwa mem -M -t 12 $As $Fr $Rr | samtools view -S -b - > "$Prefix".bam

### Add group and sample name (Prefix)
bamaddrg -b "$Prefix".bam -s $Prefix -r $Prefix > "$Prefix"_unsorted.bam
### Sort the full BAM file.
samtools sort -@ 16  "$Prefix"_unsorted.bam -o "$Prefix"_sorted.bam

#index
samtools index "$Prefix"_sorted.bam

mkdir -p $CWD/$OutDir
cp -r *_sorted* $CWD/$OutDir/.
rm -r $WorkDir
