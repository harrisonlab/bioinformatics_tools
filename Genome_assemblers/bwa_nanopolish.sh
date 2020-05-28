#!/usr/bin/env bash
#SBATCH -J canu
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=40

# Alignment of minion reads to a minion assembly prior to running nanopolish variants

Usage="sub_bwa_nanopolish.sh <assembly.fa> <reads.fa.gz> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

AssemblyIn=$1
ReadsIn=$2
OutDir=$3

echo "Assembly - $AssemblyIn"
echo "Fasta reads - $ReadsIn"
echo "OutDir - $OutDir"

CurDir=$PWD
WorkDir=$TMPDIR/nanopolish
mkdir -p $WorkDir
cd $WorkDir

Assembly=$(basename $AssemblyIn)
Reads=$(basename $ReadsIn)
cp $CurDir/$AssemblyIn $Assembly
cp $CurDir/$ReadsIn $Reads

mkdir -p $CurDir/$OutDir


bwa index $Assembly
bwa mem -x ont2d -t 16 $Assembly $Reads > aligned.sam
samtools view  -@ 8 -b aligned.sam | samtools sort  -@ 8 - -o reads.sorted.bam
# bwa mem -x ont2d -t 16 $Assembly $Reads | samtools sort -@ 16 -o reads.sorted.bam -T reads.tmp -
samtools index reads.sorted.bam

# nanopolish eventalign \
#   --progress \
#   --reads $Reads \
#   --bam reads.sorted.bam \
#   --genome $Assembly \
#   -t 16 --sam \
#   | samtools view -Su - \
#   | samtools sort - -o tmp \
#   > reads.eventalign.sorted.bam

rm $Assembly
rm $Reads
# rm aligned.sam

cp -r $WorkDir/* $CurDir/$OutDir/.
