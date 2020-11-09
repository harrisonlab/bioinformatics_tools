#!/usr/bin/env bash
#SBATCH -J bowtie
#SBATCH --partition=short
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=40


# Align raw reads to an assembly.

# ---------------
# Step 1
# Collect inputs
# ---------------

Assembly=$(basename $1)
Reads=$(basename $2)
OutDir=$3
CurDir=$PWD

if [ $4 ]; then
  AddOptions="--rg-id $4"
else
  AddOptions=""
fi


echo  "Running Bowtie with the following inputs:"
echo "Assembly - $Assembly"
echo "Trimmed reads - $Reads"
echo "OutDir - $OutDir"

# ---------------
# Step 2
# Copy data
# ---------------

WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir
cp $CurDir/$1 $Assembly
cp $CurDir/$2 $Reads

# ---------------
# Step 3
# Align seq reads
# ---------------
# Prepare the assembly for alignment
# Align reads against the assembly
# Convert the SAM file to BAM in preparation for sorting.
# Sort the BAM file, in preparation for SNP calling:
# Index the bam file

bowtie2-build $Assembly $Assembly.indexed
bowtie2 -p 16 $AddOptions -x $Assembly.indexed -U $Reads -S "$Assembly"_aligned.sam 2>&1 | tee bowtie_log.txt
samtools view --threads 8 -bS "$Assembly"_aligned.sam -o "$Assembly"_aligned.bam
samtools sort --threads 8 -o "$Assembly"_aligned_sorted.bam "$Assembly"_aligned.bam
samtools index -@ 8 "$Assembly"_aligned_sorted.bam "$Assembly"_aligned_sorted.bam.index


# ---------------
# Step 4
# Determine RPKM
# ---------------
# Determine the number of reads aligning per kilobase,
# normalised by the number of million reads aligned to
# the genome.

samtools view  --threads 8 -f 2 "$Assembly"_aligned_sorted.bam | cut -f3 | uniq -c > "$Assembly"_RPK.txt

# ---------------
# Step 5
# Cleanup
# ---------------
# Delete uneccessary files
# and copy to $OutDir

rm "$Assembly"_aligned.bam
rm *.sam
rm $Assembly
rm $Reads
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.
