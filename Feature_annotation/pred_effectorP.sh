#!/usr/bin/env bash
#SBATCH -J effectorP
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8

CurPath=$PWD
InFile="$1"
BaseName="$2"
OutDir=$CurPath/"$3"
version=$4

Organism=$(echo $InFile | rev | cut -d "/" -f3 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f2 | rev)
InName=$(echo $InFile | rev | cut -d "/" -f1 | rev)

WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile proteins.fa
echo "Running effectorP for: $Organism - $Strain"

if [ $version == "v2" ]; then
   EffectorP.py -o "$BaseName".txt -E "$BaseName".fa -i proteins.fa
  elif [ $type == "v3" ]; then
    /scratch/software/EffectorP-3.0/EffectorP-3.0-main/EffectorP.py -f -o "$BaseName".txt -E "$BaseName".fa -i proteins.fa
  else
    echo "Version 2.0 or 3.0 required"
  fi

rm proteins.fa
mkdir -p $OutDir
cp $WorkDir/* $OutDir/.
rm -r $WorkDir
