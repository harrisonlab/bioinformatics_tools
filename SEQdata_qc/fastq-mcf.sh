#!/usr/bin/env bash
#SBATCH -J fastq-mcf
#SBATCH --partition=short
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1

# Script to prepare rna for downstream applications.
# Will filter poor quality reads, perform trimming and
# remove illumina adapters.
# To be run from the project directory. Usage:
# rna_qc_fastq-mcf <RNASeq_F.fq> <RNASeq_R.fq> <illumina_adapters.fa> [DNA/RNA]

#######  Step 1	 ########
# Initialise values	#
#########################


CUR_PATH=$PWD
WORK_DIR=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

F_IN=$CUR_PATH/$1
R_IN=$CUR_PATH/$2
ILLUMINA_ADAPTERS=$3
SEQ_TYPE=$4
OutDir=$5


LIBRARY_TYPE=$(echo $F_IN | rev | cut -d "/" -f5 | rev)
ORGANISM=$(echo $F_IN | rev | cut -d "/" -f4 | rev)
STRAIN=$(echo $F_IN | rev | cut -d "/" -f3 | rev)

F_FILE=$(echo $F_IN | rev | cut -d "/" -f1 | rev | sed 's/.gz//')
R_FILE=$(echo $R_IN | rev | cut -d "/" -f1 | rev | sed 's/.gz//')

F_OUT=$(echo "$F_FILE" | sed 's/.fq/_trim.fq/g' | sed 's/.fastq/_trim.fq/g')
R_OUT=$(echo "$R_FILE" | sed 's/.fq/_trim.fq/g' | sed 's/.fastq/_trim.fq/g')

if [ $5 ]; then
  OUTDIR=$CUR_PATH/$5
else
	SEQ_TYPE=$(echo "$SEQ_TYPE" | tr "[:upper:]" "[:lower:]")
	OUTDIR="$CUR_PATH"/qc_"$SEQ_TYPE"/"$LIBRARY_TYPE"/"$ORGANISM"/"$STRAIN"
fi

echo "your compressed forward read is: $F_IN"
echo "your compressed reverse read is: $R_IN"
	echo ""
echo "your forward read is: $F_FILE"
echo "your reverse read is: $R_FILE"
	echo ""
echo "illumina adapters are stored in the file: $ILLUMINA_ADAPTERS"
	echo ""
echo "you are providing Sequence data as (DNA/RNA): $SEQ_TYPE"


#######  Step 2	 ########
# 	unzip reads			#
#########################

mkdir -p "$WORK_DIR"/F "$WORK_DIR"/R
cd "$WORK_DIR"

cat "$F_IN" | gunzip -fc > "$F_FILE"
cat "$R_IN" | gunzip -fc > "$R_FILE"

#######  Step 4	 ########
# 	Quality trim		#
#########################

fastq-mcf $ILLUMINA_ADAPTERS $F_FILE $R_FILE -o F/"$F_OUT" -o R/"$R_OUT" -C 1000000 -u -k 20 -t 0.01 -q 30 -p 5

#gzip "$WORK_DIR/*/$QC_OUTFILE"_*

gzip F/"$F_OUT"
gzip R/"$R_OUT"
mkdir -p "$OUTDIR"/F
mkdir -p "$OUTDIR"/R
cp -r F/"$F_OUT".gz "$OUTDIR"/F/"$F_OUT".gz
cp -r R/"$R_OUT".gz "$OUTDIR"/R/"$R_OUT".gz

#cat F/"$QC_OUTFILE"_F.fastq | gzip -cf > $CUR_PATH/qc_$SEQ_TYPE/paired/$ORGANISM/$STRAIN/F/"$QC_OUTFILE"_F.fastq.gz
#cat R/"$QC_OUTFILE"_R.fastq | gzip -cf > $CUR_PATH/qc_$SEQ_TYPE/paired/$ORGANISM/$STRAIN/R/"$QC_OUTFILE"_R.fastq.gz
#cp -r $WORK_DIR/F/* $CUR_PATH/qc_$SEQ_TYPE/paired/$ORGANISM/$STRAIN/F/.
#cp -r $WORK_DIR/R/* $CUR_PATH/qc_$SEQ_TYPE/paired/$ORGANISM/$STRAIN/R/.

#######  Step 8  ########
#       Cleanup         #
#########################

rm -r $WORK_DIR/
