#!/usr/bin/env bash
#SBATCH -J bowtie
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

# Align raw reads to an assembly.

# ---------------
# Step 1
# Collect inputs
# ---------------

Assembly=$(basename $1)
Read_F=$(basename $2)
Read_R=$(basename $3)
OutDir=$4

if [ $5 ]; then
  AddOptions="--rg-id $5"
else
  AddOptions=""
fi

CurDir=$PWD
echo  "Running Bowtie with the following inputs:"
echo "Assembly - $Assembly"
echo "Forward trimmed reads - $Read_F"
echo "Reverse trimmed reads - $Read_R"
echo "OutDir - $OutDir"

# ---------------
# Step 2
# Copy data
# ---------------

WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir
cp $CurDir/$1 $Assembly
cp $CurDir/$2 $Read_F
cp $CurDir/$3 $Read_R


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
bowtie2 -p 8 -X 1200 --no-mixed $AddOptions -x $Assembly.indexed  -1 $Read_F -2 $Read_R  -S "$Assembly"_aligned.sam 2>&1 | tee bowtie_log.txt

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

samtools view --threads 8 -f 2 -F 104 "$Assembly"_aligned_sorted.bam | cut -f3 | uniq -c > "$Assembly"_RPK.txt

# ---------------
# Step 5
# Cleanup
# ---------------
# Delete uneccessary files
# and copy to $OutDir

rm $Assembly
rm $Read_F
rm $Read_R
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.
rm -r $WorkDir
