#!/bin/bash
# split4signalP.sh

CurPath=$PWD
InFile=$1
Organism=$(echo $InFile | rev | cut -d "/" -f4 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f3 | rev)
InName=$(basename $InFile)

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation

SplitDir=gene_pred/final_genes_split/$Organism/$Strain
mkdir -p $SplitDir
$ProgDir/splitfile_500.py --inp_fasta $InFile --out_dir $SplitDir --out_base "$InName"_split

for File in $(ls $SplitDir/final_genes_combined.pep.fasta_split_10*); do
sbatch $ProgDir/pred_signalP.sh $File signalp-5.0
done

