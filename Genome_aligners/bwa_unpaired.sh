#!/usr/bin/env bash
#SBATCH -J bwa
#SBATCH --partition=long
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16

# Align PacBio reads to an assembly.

# ---------------
# Step 1
# Collect inputs
# ---------------

Assembly=$(basename $1)
Reads=$(basename $2)
OutDir=$3
type=$4
CurDir=$PWD
echo  "Running BWA-mem with the following inputs:"
echo "Assembly - $Assembly"
echo "Reads are - $type" 
echo "ReadsDir - $Reads"
echo "OutDir - $OutDir"


# ---------------
# Step 2
# Copy data
# ---------------

WorkDir=$TMPDIR/bwa
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

bwa index $Assembly

if [ $type == "pacbio" ]; then
bwa mem -t 24 -x pacbio $Assembly $Reads > "$Assembly"_aligned.sam
elif [ $type == "ont" ]; then
bwa mem -t 24 -x ont2d $Assembly $Reads > "$Assembly"_aligned.sam
else
echo "data type not supported"
fi

samtools view -bS "$Assembly"_aligned.sam -o "$Assembly"_aligned.bam
samtools sort "$Assembly"_aligned.bam -o "$Assembly"_aligned_sorted.bam
samtools index "$Assembly"_aligned_sorted.bam

# ---------------
# Step 4
# Show coverage
# ---------------

AlignedBam="$Assembly"_aligned_sorted.bam
CoverageTxt="$Assembly"_bp_genome_cov.txt
bedtools genomecov -max 5 -d -ibam $AlignedBam -g $Assembly > $CoverageTxt

# ---------------
# Step 5
# Flag low coverage regions
# ---------------

Threshold=5
FlaggedRegions="$Assembly"_flagged_regions.txt
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
$ProgDir/flag_low_coverage.py --genomecov $CoverageTxt --min $Threshold > $FlaggedRegions


# ---------------
# Step 6
# Cleanup
# ---------------
# Delete uneccessary files
# and copy to $OutDir

rm $Assembly
rm $Reads
rm "$Assembly"_aligned.sam
rm $CoverageTxt
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.
