#!/usr/bin/env bash
#SBATCH -J blast_pipe
#SBATCH --partition=long
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4

# script to run blast homology pipe
USAGE="blast_pipe.sh <query.fa> <dna, protein (query_format)> <genome_sequence.fa> <output_directory>"


#-------------------------------------------------------
# 		Step 0.		Initialise values
#-------------------------------------------------------

IN_QUERY=$1
QUERY_FORMAT=$2
IN_GENOME=$3
ORGANISM=$(echo $IN_GENOME | rev | cut -d "/" -f4 | rev)
STRAIN=$(echo $IN_GENOME | rev | cut -d "/" -f3 | rev)
QUERY=$(echo $IN_QUERY | rev | cut -d "/" -f1 | rev)
GENOME=$(echo $IN_GENOME | rev | cut -d "/" -f1 | rev)
CUR_PATH=$PWD
if [ "$4" ]; then OutDir=$CUR_PATH/$4; else OutDir=$CUR_PATH/analysis/blast_homology/$ORGANISM/$STRAIN; fi
WORK_DIR=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir $WORK_DIR
cd $WORK_DIR
cp $CUR_PATH/$IN_GENOME $GENOME
cp $CUR_PATH/$IN_QUERY $QUERY
OUTNAME="$STRAIN"_"$QUERY"
SCRIPT_DIR=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis

echo "Running blast_pipe.sh"
echo "Usage = $USAGE"
echo "Organism is: $ORGANISM"
echo "Strain is: $STRAIN"
echo "Query is: $QUERY"
echo "This is $QUERY_FORMAT data"
echo "Genome is: $GENOME"
echo "You are running scripts from:"
echo "$SCRIPT_DIR"

if test "$QUERY_FORMAT" = 'protein'; then
	SELF_BLAST_TYPE='blastp'
	BLAST_CSV_TYPE='tblastn'
elif test "$QUERY_FORMAT" = 'dna'; then
	SELF_BLAST_TYPE='blastn'
	BLAST_CSV_TYPE='tblastx'
#	BLAST_CSV_TYPE='blastn'
else exit
fi


#-------------------------------------------------------
# 		Step 1.		blast queries against themselves
#-------------------------------------------------------

$SCRIPT_DIR/blast_self.pl $QUERY $SELF_BLAST_TYPE > "$QUERY"_self.csv

#-------------------------------------------------------
# 		Step 2.		simplify hits table into homolog groups
#-------------------------------------------------------

$SCRIPT_DIR/blast_parse.pl "$QUERY"_self.csv > "$QUERY"_simplified.csv

#-------------------------------------------------------
# 		Step 3.		blast queries against genome
#-------------------------------------------------------

$SCRIPT_DIR/blast2csv.pl $QUERY $BLAST_CSV_TYPE $GENOME 5 > "$OUTNAME"_hits.csv

#-------------------------------------------------------
# 		Step 4.		combine the homolog group table
#					 with the blast result table
#-------------------------------------------------------

paste -d '\t' "$QUERY"_simplified.csv <(cut -f 2- "$OUTNAME"_hits.csv) > "$OUTNAME"_homologs.csv

#-------------------------------------------------------
# 		Step 5.		Cleanup
#-------------------------------------------------------

mkdir -p $OutDir

cp -r $WORK_DIR/"$OUTNAME"_homologs.csv $OutDir/.

rm -r $WORK_DIR/
