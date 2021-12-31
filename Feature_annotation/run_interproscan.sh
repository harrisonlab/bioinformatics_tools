#!/bin/bash
#SBATCH -J interproscan
#SBATCH --partition=long
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=10

# run_interproscan.sh
# The current version of interproscan only works with Java version 11
# Note - the latest version of interproscan doesnt work on the cluster
# because of GLIBC cant be updated.
# PATH=/data/scratch/gomeza/prog/java/jdk-11.0.4/bin:${PATH}

CUR_PATH=$PWD
IN_FILE=$1
ORGANISM=$(echo $IN_FILE | rev | cut -d "/" -f3 | rev)
STRAIN=$(echo $IN_FILE | rev | cut -d "/" -f2 | rev)
IN_NAME=$(basename $IN_FILE)

WORK_DIR=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WORK_DIR
cd $WORK_DIR
cp $CUR_PATH/$IN_FILE $IN_NAME
sed -i -r 's/\*/X/g' $IN_NAME

/scratch/software/interproscan-5.52-86.0/interproscan.sh -appl CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,SFLD,SMART,SUPERFAMILY,TIGRFAM -goterms -iprlookup -pa -i $IN_NAME

OUT_DIR=$CUR_PATH/gene_pred/interproscan/$ORGANISM/$STRAIN/raw
mkdir -p $OUT_DIR
cp *.gff3 $OUT_DIR/.
cp *.tsv $OUT_DIR/.
cp *.xml $OUT_DIR/.

#rm -r $WORK_DIR