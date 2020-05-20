#!/bin/bash
# sub_interproscan.sh

CurPath=$PWD
InFile=$1
Organism=$(echo $InFile | rev | cut -d "/" -f4 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f3 | rev)
InName=$(basename $InFile)

InterproDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides

SplitDir=gene_pred/interproscan/$Organism/$Strain
mkdir -p $SplitDir
$SplitfileDir/splitfile_500.py --inp_fasta $InFile --out_dir $SplitDir --out_base "$InName"_split
for file in $(ls $SplitDir/*_split*); do
	Jobs=$(qstat | grep 'run_interp' | grep 'qw' | wc -l)
	while [ $Jobs -gt 1 ]; do
		sleep 10
		printf "."
		Jobs=$(qstat | grep 'run_interp' | grep 'qw' | wc -l)
	done
	printf "\n"
	echo $file
	qsub $InterproDir/run_interproscan.sh $file
done
