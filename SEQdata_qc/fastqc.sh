#!/usr/bin/env bash
#SBATCH -J fastqc
#SBATCH --partition=short
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1

INFILE=$1
DATA_TYPE=$(echo $INFILE | rev | cut -d "/" -f6 | rev)
READ_TYPE=$(echo $INFILE | rev | cut -d "/" -f5 | rev)
ORGANISM=$(echo $INFILE | rev | cut -d "/" -f4 | rev)
STRAIN=$(echo $INFILE | rev | cut -d "/" -f3 | rev)
DIRECTION=$(echo $INFILE | rev | cut -d "/" -f2 | rev)
READS=$(echo $INFILE | rev | cut -d "/" -f1 | rev)

CUR_PATH=$PWD
WORK_DIR=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WORK_DIR
cd $WORK_DIR

cp $CUR_PATH/$INFILE .

fastqc --nogroup $READS

unzip *.zip

#mkdir -p $CUR_PATH/$DATA_TYPE/$READ_TYPE/$ORGANISM/$STRAIN/$DIRECTION/

cp -r $WORK_DIR/*fastqc $CUR_PATH/$DATA_TYPE/$READ_TYPE/$ORGANISM/$STRAIN/$DIRECTION/.

rm -r $WORK_DIR

exit