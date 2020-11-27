#!/usr/bin/env bash
#SBATCH -J kat
#SBATCH --partition=short
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=24

# Assess assembly quality by plotting occurence of kmers in a genome assembly
# vs their occurence in illumina data.


Threads=8

# ---------------
# Step 1
# Collect inputs
# ---------------

Assembly=$(basename $1)
Read_F=$(basename $2)
Read_R=$(basename $3)
OutDir=$4
Prefix=$5
if [ $6 ]; then
  MaxX=$6
fi


CurDir=$PWD
echo "Running KAT with the following inputs:"
echo "Assembly - $Assembly"
echo "Forward reads - $Read_F"
echo "Reverse reads - $Read_R"
echo "OutDir - $OutDir"
echo "Outfile prefix = $Prefix"

# ---------------
# Step 2
# Copy data
# ---------------

WorkDir=$TMPDIR/kat
mkdir -p $WorkDir
cd $WorkDir
cp $CurDir/$1 $Assembly
# cp $CurDir/$2 $Read_F
# cp $CurDir/$3 $Read_R
cat $CurDir/$2 $CurDir/$3 > reads.fa

# ---------------
# Step 3
# Generate KAT kmer spectrum
# ---------------

mkdir out
kat comp -o out/${Prefix} -m 21 -t 8 reads.fa $Assembly

ProgDir=/home/gomeza/miniconda3/envs/qc_tools/share/kat/scripts/kat/plot
$ProgDir/spectra_cn.py -o out/${Prefix}_spectra.png out/${Prefix}-main.mx
if [ $MaxX ]; then
  $ProgDir/spectra_cn.py -o out/${Prefix}_spectra_X-${MaxX}.png -x $MaxX out/${Prefix}-main.mx
fi

# ---------------
# Step 4
# Kmer coverage per contig
# ---------------

kat sect -o out/${Prefix}_sect_reads -m 21 -t 8 $Assembly reads.fa
kat sect -o out/${Prefix}_sect_assembly -m 21 -t 8 $Assembly $Assembly


mkdir -p ${CurDir}/${OutDir}
mv out/* ${CurDir}/${OutDir}/.
