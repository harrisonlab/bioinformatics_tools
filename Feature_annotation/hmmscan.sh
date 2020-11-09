#!/usr/bin/env bash
#SBATCH -J tmhmm
#SBATCH --partition=long
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8

# hmmscan.sh <HmmFile> <ProtFile> <OutputPrefix> <OutputDirecotry>

CurPath=$PWD
HmmFile=$1
ProtFile=$2
Prefix=$3
OutDir=$4


WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WorkDir
cd $WorkDir
cp $CurPath/$ProtFile proteins.fa

hmmscan --cpu 8 --domtblout $Prefix.out.dm $CurPath/$HmmFile proteins.fa > $Prefix.out

mkdir -p $CurPath/$OutDir
cp $Prefix.out.dm $CurPath/$OutDir/.
cp $Prefix.out $CurPath/$OutDir/.

rm -r $WorkDir