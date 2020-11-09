#!/bin/bash
# append_interpro.sh

USAGE="append_interpro.sh predicted_genes.fa path_to_interpro_results"

CUR_PATH=$PWD
GENE_MODELS=$1
RESULTS_PATH=$2
ORGANISM=$(echo $GENE_MODELS | rev | cut -d "/" -f4 | rev)
STRAIN=$(echo $GENE_MODELS | rev | cut -d "/" -f3 | rev)
IN_NAME=$(basename "$IN_FILE")

OUT_PATH=gene_pred/interproscan/$ORGANISM/$STRAIN

mkdir -p $OUT_PATH
printf "" > $OUT_PATH/"$STRAIN"_interproscan.tsv
printf "" > $OUT_PATH/"$STRAIN"_interproscan.xml
printf "" > $OUT_PATH/interpro_features.gff
printf "" > $OUT_PATH/"$STRAIN"_interpro.gff3
for FILE in $(ls -v $RESULTS_PATH/*_split_*.tsv); do
	cat $FILE >> $OUT_PATH/"$STRAIN"_interproscan.tsv
done
for FILE in $(ls -v $RESULTS_PATH/*_split_*.xml); do
	cat $FILE >> $OUT_PATH/"$STRAIN"_interproscan.xml
done
for FILE in $(ls -v $RESULTS_PATH/*_split_*.gff3); do
	FILE_NAME=$(basename $FILE)
	FASTA_START=$(cat $FILE | grep -E "^##FASTA" -n | cut -d ':' -f1)
	cat $FILE | head -n "$FASTA_START" | grep -v -E "^#" >> $OUT_PATH/interpro_features.gff
done
cat $OUT_PATH/interpro_features.gff $GENE_MODELS >> $OUT_PATH/"$STRAIN"_interpro.gff3
rm $OUT_PATH/interpro_features.gff
