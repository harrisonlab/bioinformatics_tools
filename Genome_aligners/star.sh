#!/usr/bin/env bash
#SBATCH -J STAR
#SBATCH --partition=long
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=16

#Align RNAseq data with genome using STAR

Usage="star.sh InGenome.fa InReadF.fa InReadR.fa Output_directory [GeneLocations.gff]"
echo $Usage

# ---------------
# Step 1
# Collect inputs
# ---------------
GffProvided="N"

InGenome=$(basename $1)
InReadF=$(basename $2)
InReadR=$(basename $3)
OutDir=$4

# determine if optional file for genemodels has been provided
if [ $5 ]; then
  GffProvided="Y"
  InGff=$(basename $5)
fi

# Set working directory
CurDir=$PWD
WorkDir=$TMPDIR/star
GenomeDir=$WorkDir/index
mkdir -p $GenomeDir


# Copy over input files
cd $WorkDir

cp $CurDir/$1 $InGenome
cp $CurDir/$2 $InReadF
cp $CurDir/$3 $InReadR

if [ $GffProvided == "Y" ]; then
cp $CurDir/$5 $InGff
fi


# ---------------
# Step 2
# Create the Index File
# ---------------
echo "Building index file"
ParentFeature="Parent"


if [ $GffProvided == "N" ]; then
STAR \
--runMode genomeGenerate \
--genomeDir $GenomeDir \
--genomeFastaFiles $InGenome \
--runThreadN 8
elif [ $GffProvided == "Y" ]; then
STAR \
--runMode genomeGenerate \
--genomeDir $GenomeDir \
--genomeFastaFiles $InGenome \
--sjdbGTFtagExonParentTranscript $ParentFeature \
--sjdbGTFfile $InGff \
--runThreadN 8 \
--sjdbOverhang 99
fi

# ---------------
# Step 2=3
# Run STAR
# ---------------

echo "Aligning RNAseq reads"

STAR \
--genomeDir $GenomeDir \
--outFileNamePrefix star_aligment \
--readFilesCommand zcat \
--readFilesIn $InReadF $InReadR \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--runThreadN 8


#STAR \
#--genomeDir $GenomeDir \
#--outFileNamePrefix star_aligment \
#--readFilesCommand zcat \
#--readFilesIn $InReadF $InReadR \
#--outSAMtype BAM SortedByCoordinate \
#--outSAMstrandField intronMotif \
#--winAnchorMultimapNmax 200 \
#--seedSearchStartLmax 30 \
#--quantMode TranscriptomeSAM GeneCounts \
#--outReadsUnmapped Fastx \
#--runThreadN 8

# --genomeDir $GenomeDir \
# --outFileNamePrefix star_aligment \
# --readFilesIn $InReadF $InReadR \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMstrandField intronMotif \
# --runThreadN 8


rm -r $GenomeDir
rm $InGenome
rm $InGff
rm $InReadF
rm $InReadR
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.
