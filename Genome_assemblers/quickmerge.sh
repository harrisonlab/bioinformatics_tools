#!/usr/bin/env bash
#SBATCH -J quickmerge
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16

# Merge a Pacbio only (Canu) and a hybrid pacbio-illumina (Spades) assembly
# to make a merged assembly.
# Submit the best assembly as the first input and the worse of the two Assemblies
# as the second input. The first input will be used to scaffold.
# The N50 is recomended to be used as the Anchor length for contigs. However,
# beware that if the N50 is very large (>1Mb) then problems may be encountered
# if trying to scaffold small chromosomes (<1Mb).

Usage="quickmerge.sh <better_assembly.fa> <worse_assembly.fa> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

PacbioFa=$1
HybridFa=$2
OutDir=$3
Prefix=pacbio_hybrid_merged
AnchorOverlapCuttoff=5.0
ExtensionOverlapCuttoff=1.5
AnchorContigLength=$4
MinAlignmentLength=5000

echo "Running Quickmerge with the following inputs:"
echo "Pacbio assembly (or best quality assembly) - $PacbioFa"
echo "Hybrid assembly (or worse assembly) - $HybridFa"
echo "OutDir - $OutDir"
echo "Contigs of this size will be used as anchors for merging - $AnchorContigLength"

# ---------------
# Step 2
# Copy data
# ---------------

CurPath=$PWD
WorkDir="$TMPDIR"/quickmerge
mkdir -p $WorkDir
cd $WorkDir
cp $CurPath/$PacbioFa pacbio.fa
cp $CurPath/$HybridFa hybrid.fa

# ---------------
# Step 3
# Run quickmerge
# ---------------

merge_wrapper.py \
  pacbio.fa \
  hybrid.fa \
  -pre pacbio_hybrid_merged \
  -hco $AnchorOverlapCuttoff \
  -c $ExtensionOverlapCuttoff \
  -l $AnchorContigLength \
  -ml $MinAlignmentLength \
  > merge_report.txt

rm pacbio.fa
rm hybrid.fa

mkdir -p $CurPath/$OutDir
# mv merged.fasta $CurPath/$OutDir/.
# mv merge_report.txt $CurPath/$OutDir/.
cp -r $WorkDir/* $CurPath/$OutDir/.

rm -r $WorkDir
