#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace01.blacklace|blacklace05.blacklace|blacklace06.blacklace


CurPath=$PWD
InFile="$1"
BaseName="$2"
OutDir=$CurPath/"$3"

Organism=$(echo $InFile | rev | cut -d "/" -f3 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f2 | rev)
InName=$(echo $InFile | rev | cut -d "/" -f1 | rev)

WorkDir=$TMPDIR/"$Strain"_"$InName"
mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile proteins.fa
echo "Running effectorP for: $Organism - $Strain"
EffectorP.py -o "$BaseName".txt -E "$BaseName".fa -i proteins.fa

rm proteins.fa
mkdir -p $OutDir
cp $WorkDir/* $OutDir/.
