#!/usr/bin/env bash
#SBATCH -J ORF
#SBATCH --partition=long
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8

set -eu
set -o pipefail

Usage='ORF_finder.sh <Assembly.fa>'


#######  Step 1	 ########
#   Initialise values	#
#########################

CurPath=$PWD
Assembly=$1

ScriptDir=/home/agomez/scratch/apps/scripts/Feature_annotation

#Organism=$(echo $Assembly | rev | cut -d "/" -f4 | rev)
#Strain=$(echo $Assembly | rev | cut -d "/" -f6 | rev)
#Organism=N.ditissima
Strain=Bra_1
SortedContigs=$(basename $Assembly)

WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WorkDir
cd $WorkDir
cp $CurPath/$Assembly $SortedContigs

echo $SortedContigs
echo "The following files are present in the temporary directory:"
ls

echo $ScriptDir


#######  Step 2 ########
# find open reading frames in contigs and translate these to AA	#
#########################

echo "Predicting coding seqs- forward"
$ScriptDir/print_atg_gff.pl $SortedContigs F "$Strain"_F_atg.fa "$Strain"_F_atg_nuc.fa "$Strain"_F_atg_ORF.gff

#######  Step 3 ########
# revcomp contigs to get reads on the R strand of DNA #
#########################

echo "REVCOMPing the contigs"
$ScriptDir/revcomp_fasta.pl $SortedContigs > contigs_R.fa


#######  Step 4 ########
# find open reading frames in revcomp contigs and translate	these to AA #
#########################

echo "Predicting coding seqs- reverse"
$ScriptDir/print_atg_gff.pl contigs_R.fa R "$Strain"_R_atg.fa "$Strain"_R_atg_nuc.fa "$Strain"_R_atg_ORF.gff

#######  Step 4 ########
# concatenate files of open reading frames #
#########################

echo "Joining Forward and Reverse Files"
cat "$Strain"_F_atg.fa "$Strain"_R_atg.fa > $Strain.aa_cat.fa
cat "$Strain"_F_atg_nuc.fa "$Strain"_R_atg_nuc.fa > "$Strain"_nuc.fa
cat "$Strain"_F_atg_ORF.gff "$Strain"_R_atg_ORF.gff > "$Strain"_ORF.gff

#######  Step 6 ########
#       Cleanup         #
#########################

echo "Cleaning up any old files"
rm -rf ./fasta_seq

mkdir -p $CurPath/gene_pred/ORF_finder/$Organism/$Strain/
cp -r $WorkDir/. $CurPath/gene_pred_vAG/ORF_finder/$Organism/$Strain/.
rm -r $TMPDIR

