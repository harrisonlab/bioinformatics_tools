#!/usr/bin/env bash
#SBATCH -J transposonPSI
#SBATCH --partition=short
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=60

Usage="sbatch transposonPSI.sh <assembled_contigs.fa> [<output_directory>]"

InFile=$1

Organism=$(echo $InFile | rev | cut -d "/" -f4 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f3 | rev)
Assembly=$(echo $InFile | rev | cut -d "/" -f5 | rev)

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

if [ $2 ]; then
  OutDir=$CurPath/$2
else
  OutDir=$CurPath/repeat_masked/$Organism/$Strain/"$Assembly"_repmask
fi

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile "$Strain"_contigs_unmasked.fa

# The full path is needed
/home/gomeza/miniconda3/envs/general_tools/share/transposonPSI/transposonPSI.pl "$Strain"_contigs_unmasked.fa nuc

mkdir -p $OutDir
cp -r $WorkDir/* $OutDir/.

rm -r $WorkDir
