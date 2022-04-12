#!/usr/bin/env bash
#SBATCH -J effectorP
#SBATCH --partition=long
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

if [ $version == "2.0" ]; then
  /home/agomez/scratch/apps/prog/EffectorP/EffectorP_2.0/Scripts/EffectorP.py -o "$BaseName".txt -E "$BaseName".fa -i proteins.fa
elif [ $version == "3.0" ]; then
  python /home/agomez/scratch/apps/prog/EffectorP/EffectorP_3.0.0-beta/EffectorP.py -f -o "$BaseName".txt -E "$BaseName".fa -i proteins.fa
else
 echo "Version 2.0 or 3.0"
fi

rm proteins.fa
mkdir -p $OutDir
cp $WorkDir/* $OutDir/.
rm -r $WorkDir
