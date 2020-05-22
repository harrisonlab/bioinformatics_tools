#!/usr/bin/env bash
#SBATCH -J STAR
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=12G
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

echo $InGenome
echo $InReadF
echo $InReadR

# determine if optional file for genemodels has been provided
if [ $5 ]; then
  GffProvided="Y"
  InGff=$(basename $5)
fi

# Set working directory
CurDir=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
echo $WorkDir
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
--genomeSAindexNbases 11 \
--genomeDir $GenomeDir \
--genomeFastaFiles $InGenome \
--runThreadN 16
elif [ $GffProvided == "Y" ]; then
STAR \
--runMode genomeGenerate \
--genomeSAindexNbases 11 \
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
--outSAMtype BAM Unsorted \
--outSAMstrandField intronMotif \
--runThreadN 16


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

samtools sort star_aligmentAligned.out.bam > star_aligmentAligned.sorted.out.bam

rm -r $GenomeDir
rm $InGenome
rm $InGff
rm $InReadF
rm $InReadR
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.

rm -r $WorkDir