#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l virtual_free=1G
##$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace

# run_interproscan.sh
# The current version of interproscan only works with Java version 11
# Note - the latest version of interproscan doesnt work on the cluster
# because of GLIBC cant be updated.
# PATH=/home/armita/prog/java/jdk-11.0.4/bin:${PATH}

CUR_PATH=$PWD
IN_FILE=$1
ORGANISM=$(echo $IN_FILE | rev | cut -d "/" -f3 | rev)
STRAIN=$(echo $IN_FILE | rev | cut -d "/" -f2 | rev)
IN_NAME=$(basename $IN_FILE)

WORK_DIR=$TMPDIR/interpro

mkdir -p $WORK_DIR
cd $WORK_DIR
cp $CUR_PATH/$IN_FILE $IN_NAME
sed -i -r 's/\*/X/g' $IN_NAME

PYTHONPATH="/home/armita/.local/lib/python3.5/site-packages:/home/armita/prog/kat"
interproscan.sh -goterms -iprlookup -pa -i $IN_NAME

OUT_DIR=$CUR_PATH/gene_pred/interproscan/$ORGANISM/$STRAIN/raw
mkdir -p $OUT_DIR
cp *.gff3 $OUT_DIR/.
cp *.tsv $OUT_DIR/.
cp *.xml $OUT_DIR/.
