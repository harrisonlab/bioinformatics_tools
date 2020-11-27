#!/usr/bin/env bash
#SBATCH -J mummer
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8


# Whole genome alignment

# ---------------
# Step 1
# Collect inputs
# ---------------

SubjectGenome=$(basename $1)
QueryGenome=$(basename $2)
Prefix=$3
OutDir=$4

CurDir=$PWD
echo  "Running MUMmer with the following inputs:"
echo "Subject genome - $SubjectGenome"
echo "Query genome - $QueryGenome"
echo "Outputs will be labelled with Prefix: $Prefix"
echo "OutDir - $OutDir"

# ---------------
# Step 2
# Copy data
# ---------------

WorkDir=$TMPDIR/MUMmer
mkdir -p $WorkDir
cd $WorkDir
cp $CurDir/$1 $SubjectGenome
cp $CurDir/$2 $QueryGenome


# ---------------
# Step 3
# Align Query to subject
# ---------------
#

/scratch/software/mummer/mummer-4.0.0rc1/nucmer -t 8 -p $Prefix $SubjectGenome $QueryGenome > $Prefix.txt

/scratch/software/mummer/mummer-4.0.0rc1/delta-filter -i 70 -l 5000 -r -q $Prefix.delta > "$Prefix"_filtered.delta

/scratch/software/mummer/mummer-4.0.0rc1/show-tiling -p ${Prefix}_tiled_pseudomolecules.fa "$Prefix"_filtered.delta > "$Prefix".tiling

/scratch/software/mummer/mummer-4.0.0rc1/mummerplot -p $Prefix "$Prefix"_filtered.delta

/scratch/software/mummer/mummer-4.0.0rc1/show-coords -c -l -b -T "$Prefix"_filtered.delta > "$Prefix"_coords.tsv

# ---------------
# Step 5
# Cleanup
# ---------------
# Delete uneccessary files
# and copy to $OutDir

rm $SubjectGenome
rm $QueryGenome
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.
