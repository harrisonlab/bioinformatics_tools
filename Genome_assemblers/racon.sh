#!/usr/bin/env bash
#SBATCH -J racon
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=16

# Alignment of minion reads to a minion assembly prior to running nanopolish variants

Usage="sub_racon.sh <assembly.fa> <corrected_reads.fq.gz> <number_of_iterations> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

AssemblyIn=$1
ReadsIn=$2
Iterations=$3
OutDir=$4

echo "Assembly - $AssemblyIn"
echo "Fasta reads - $ReadsIn"
echo "OutDir - $OutDir"

CurDir=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

Assembly=$(basename $AssemblyIn)
Assembly=$(echo $Assembly | sed 's/.utg/.fa/g')
Reads=$(basename $ReadsIn)
cp $CurDir/$AssemblyIn $Assembly
cp $CurDir/$ReadsIn $Reads

Prefix=$(echo $Assembly | cut -f1 -d '.')

mkdir -p $CurDir/$OutDir

cp $Assembly current-assembly.fa
for i in $(seq 1 $Iterations); do
echo "Iteration - $i"
minimap2 \
-x map-ont \
-t16 \
current-assembly.fa \
$Reads \
> racon_round_$i.reads_mapped.paf

racon $Reads racon_round_$i.reads_mapped.paf current-assembly.fa > $WorkDir/racon_round_$i.fasta
cp racon_round_$i.fasta current-assembly.fa
cp racon_round_$i.fasta $CurDir/$OutDir/"$Prefix"_racon_round_$i.fasta
done

rm -r $WorkDir
