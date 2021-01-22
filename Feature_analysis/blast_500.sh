#!/usr/bin/env bash
#SBATCH -J blastall
#SBATCH --partition=long
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=8

set -u
set -e

CurPath=$PWD
BlastDB=$1
Infile=$2
OutFile=$3

echo "Infiles:"
echo $BlastDB
echo $Infile
echo $OutFile

Db=$(basename $1)
Query=$(basename $2)
Hits=$(basename $3)

echo "Copying to Files:"
echo $Db
echo $Query

echo "Outfiles will be named:"
echo $Hits

WorkDir=$TMPDIR/blastall
mkdir -p $WorkDir
cd $WorkDir

for File in $(ls $CurPath/"$BlastDB"*); do
  cp $File .
done
cp $CurPath/$Infile $Query

echo "hits will be named:"
echo $Hits
echo "hits will be moved to:"
echo "$CurPath/$OutFile"

blastall -d $Db -p blastp -i $Query -v 100000 -b 100000 -e 1e-5 -m 8 -F 'm S' -a 1 -o $Hits

cp $Hits $CurPath/$OutFile
