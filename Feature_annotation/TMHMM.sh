#!/usr/bin/env bash
#SBATCH -J tmhmm
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=8

USAGE='TMHMM.sh <predicted_proteins.fa>'
INFILE=$1

echo $USAGE

ORGANISM=$(echo $INFILE | rev | cut -d "/" -f4 | rev)
STRAIN=$(echo $INFILE | rev | cut -d "/" -f3 | rev)
PRED_PROTEINS=$(echo $INFILE | rev | cut -d "/" -f1 | rev)
CUR_PATH=$PWD
WORK_DIR=$CUR_PATH/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WORK_DIR
cd $WORK_DIR

cat $CUR_PATH/$INFILE | /home/agomez/scratch/apps/prog/tmhmm-2.0c/bin/tmhmm > "$STRAIN"_tmhmm_out.txt
cat "$STRAIN"_tmhmm_out.txt | grep -v -w 'PredHel=0' > "$STRAIN"_TM_genes_pos.txt
cat "$STRAIN"_tmhmm_out.txt | grep -w 'PredHel=0' > "$STRAIN"_TM_genes_neg.txt

mkdir -p $CUR_PATH/gene_pred/trans_mem/$ORGANISM/$STRAIN/
cp -r $WORK_DIR/* $CUR_PATH/gene_pred/trans_mem/$ORGANISM/$STRAIN/.

rm -r $WORK_DIR
