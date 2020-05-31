#!/usr/bin/env bash
#SBATCH -J nanopolish_variants
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=40

# nanopolish variants: detect SNPs and indels with respect to a reference genome
# nanopolish variants --consensus: calculate an improved consensus sequence for a draft genome assembly

Usage="nanopolish_variants.sh <assembly.fa> <reads.fa.gz> <aligned_ONTreads.sam> <ploidy[e.g.1]> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

AssemblyIn=$1
ReadsIn=$2
AlignedIn=$3
Ploidy=$4
#Region=$5
OutDir=$5
#Prefix=$(echo $Region | sed -e "s/ /_/g")

echo "Assembly - $AssemblyIn"
echo "Fasta reads - $ReadsIn"
echo "Aligned reads - $AlignedIn"
echo "OutDir - $OutDir"

CurDir=$PWD
WorkDir=$CurDir/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WorkDir
cd $WorkDir

Assembly=$(basename $AssemblyIn)
Reads=$(basename $ReadsIn)
Aligned=$(basename $AlignedIn)
cp $CurDir/$AssemblyIn $Assembly
cp $CurDir/$AlignedIn $Aligned
cp $CurDir/$AlignedIn.bai $Aligned.bai
cp $CurDir/$ReadsIn $Reads
cp $CurDir/$ReadsIn* .

echo "Files in the folder are:"
ls

nanopolish_makerange.py $Assembly | parallel --results nanopolish.results -P 8 \
nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r $Reads \
--ploidy $Ploidy \
--max-haplotypes 100000 \
--fix-homopolymers \
--min-candidate-frequency 0.2 \
-b $Aligned \
-g $Assembly -t 4

#nanopolish variants \
  #-t 8 \
  #--ploidy $Ploidy \
  #-w $Region \
  #--consensus \
  #--max-haplotypes 100000 \
  #--fix-homopolymers \
  #--min-candidate-frequency 0.2 \
  #--reads $Reads \
  #--bam $Aligned \
  #--genome $Assembly
  #> "$Prefix"_variants.txt

mkdir -p $CurDir/$OutDir
#cp "$Prefix"* $CurDir/$OutDir/.

cp -r $WorkDir/* $CurDir/$OutDir/.
rm -r $WorkDir
