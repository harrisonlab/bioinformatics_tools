#!/usr/bin/env bash
#SBATCH -J featureCounts
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

# Counting reads to genomic features such as genes, exons, promoters and genomic bins.
Usage='squeue $ProgDir/featureCounts.sh $BamFile $Gff $OutDir $OutFilePrefix'

# ---------------
# Step 1
# Collect inputs
# ---------------

InBam=$(basename $1)
InGff=$(basename $2)
OutDir=$3
Prefix=$4

CurDir=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
echo "$WorkDir"
mkdir -p $WorkDir
cd $WorkDir

cp $CurDir/$1 $InBam
cp $CurDir/$2 $InGff


# ---------------
# Step 2
# Run featureCounts
# ---------------

featureCounts \
  -p -B -R BAM \
  -M --fraction \
  -T 4 \
  -a $InGff \
  -t exon \
  -g "Parent" \
  -o "$Prefix"_featurecounts.txt \
  $InBam

rm $InBam
rm $InGff
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.

rm -r $WorkDir