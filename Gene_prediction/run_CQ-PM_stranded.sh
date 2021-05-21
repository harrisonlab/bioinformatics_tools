#!/bin/bash

#This is the prediction from transcript run script
#for stranded RNA-seq. This script automates
#CodingQuarry's pathogen mode (see manual).

#Provide a genome and gff of transcripts to run:
#run_CQ-PM_stranded.sh myGenome.fa myTranscripts.gff
#See the manual for input format requirements.

#This version requires signalp 5

echo "Running the standard CodingQuarry predictions..."
CodingQuarry -f $1 -t $2 -p 4

echo "Translating and quality filtering the coding sequence..."
python $QUARRY_PATH/scripts/fastaTranslate.py out/Predicted_CDS.fa | sed 's/*$//g' > CQ_Proteins.fa
python $QUARRY_PATH/scripts/gene_errors_Xs.py CQ_Proteins.fa CQPMtmp.fa
mv CQPMtmp.fa CQ_Proteins.fa

echo "Running signalP..."
#Split the protein file up into smaller chunks
python $QUARRY_PATH/scripts/split_fasta.py CQ_Proteins.fa 200

#Add ".fasta" to every file. 
    for file in CQ_Proteins.fa-*
    do
        mv "$file" "$file".fasta
    done

i=0
    for FILE in CQ_Proteins.fa-*
    do
        #If signalP is not in your path, change the line below to specify its location
        signalp-5.0 -fasta $FILE > CQ_Proteins_out_$i
        rm $FILE
        i=$(($i+1))
    done
cat CQ_Proteins*_summary.signalp5 | grep -v "#" | awk '($2 == "SP(Sec/SPI)"){print $1" "$7}' > Secretome.txt
rm CQ_Proteins_out_*
rm CQ_Proteins.fa

echo "Running CodingQuarry-PM..."
CodingQuarry -f $1 -t $2 -2 out/PredictedPass.gff3 -p 4 -g Secretome.txt -h



