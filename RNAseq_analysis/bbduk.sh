#!/usr/bin/env bash
#SBATCH -J bbduk
#SBATCH --partition=long
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=12
#SBATCH -o bbduk_"%j".out

# Script for bbduk program: Decontamination of rRNA reads in RNAseq data

RIBOKMERS=$1; shift
OUTDIR=$1; shift
FORWARD=$1; shift
REVERSE=$1; shift
PROGDIR=$1; shift
STRAIN=$1; shift

echo "STRAIN" $STRAIN

CUR_PATH=$PWD
echo "CUR_PATH" $CUR_PATH

#RNA=${1:-true}; shift
#TRIML=( ktrim=l k=23 mink=11 hdist=1 tpe tbo t=10 )
#TRIMR=( ktrim=r k=23 mink=11 hdist=1 tpe tbo t=10 )
#PHIX=( k=31 hdist=1 t=10 )
RRNA=( k=31 t=10 )

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

echo "tempDIR" $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

echo "FOR" $FORWARD

F=$(sed 's/.*\///' <<<$FORWARD)
R=$(sed 's/.*\///' <<<$REVERSE)

echo "F" $F
echo "R" $R

# Input files

F_IN=$CUR_PATH/$FORWARD
R_IN=$CUR_PATH/$REVERSE

echo "F_IN" $F_IN
echo "R_IN" $R_IN

# Output files

FOUT=$(sed 's/trim.fq.gz//' <<<$F)
ROUT=$(sed 's/trim.fq.gz//' <<<$R)

OUT=F/"$FOUT"cleaned.fq.gz
OUT2=R/"$ROUT"cleaned.fq.gz
OUTM=F/"$FOUT"rRNA.fq.gz
OUTM2=R/"$ROUT"rRNA.fq.gz
STATS=$STRAIN.stats.txt

echo "OUT" $OUT
echo "OUT2" $OUT2
echo "OUTM" $OUTM
echo "OUTM2" $OUTM2
echo "STATS" $STATS


bbduk.sh threads=10 in=$F_IN in2=$R_IN out=$OUT out2=$OUT2 outm=$OUTM outm2=$OUTM2 ref=$RIBOKMERS ${RRNA[@]} stats=$STATS

echo "Finish and clean"

# Move output files
mkdir -p "$CUR_PATH"/"$OUTDIR"/stats
mkdir -p "$CUR_PATH"/"$OUTDIR"/F
mkdir -p "$CUR_PATH"/"$OUTDIR"/R
cp -r F/*.fq.gz "$CUR_PATH"/"$OUTDIR"/F/.
cp -r R/*.fq.gz "$CUR_PATH"/"$OUTDIR"/R/.

echo "outdir/stats" $CUR_PATH/$OUTDIR/stats

# cp *.fq.gz $OUTDIR/.
cp $STATS "$CUR_PATH"/"$OUTDIR"/stats/.
echo "stats" $STATS

#  cleanup
cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
