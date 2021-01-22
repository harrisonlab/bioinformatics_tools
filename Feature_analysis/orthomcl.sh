#!/usr/bin/env bash
#SBATCH -J orthomcl
#SBATCH --partition=long
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=8

set -u
set -e
set -o pipefail

Usage='orthomcl.sh <merged_file_of_blast_hits.tsv> <good_proteins.fasta> [Inflation value (1-5)]'

# ----------------------	Step 1	----------------------
# 		Set Variables
#
#-------------------------------------------------------

MergeHits=$1
GoodProts=$2
# The inflation value determines the tightness of your clusters
#   - If higher you will have smaller and more closely-related clusters
Inflation='1.5'
if [ -n "$3" ]; then
  Inflation="$3"
fi
IsolateAbrv=$(echo $GoodProts | rev | cut -f3 -d '/' | rev)

CurPath=$PWD
WorkDir=$TMPDIR/orthomcl
OutDir=$CurPath/analysis/orthology/orthomcl/$IsolateAbrv
mkdir -p $WorkDir
cd $WorkDir

# HostIP='149.155.34.104'
HostIP='192.168.1.100'

Config="$IsolateAbrv"_orthomcl.config
SimilarGenes="$IsolateAbrv"_similar.txt
Log_file="$IsolateAbrv"_orthoMCL.log
MclInput="$IsolateAbrv"_mclInput
MclOutput="$IsolateAbrv"_mclOutput
OrthoGroups="$IsolateAbrv"_orthogroups.txt
OrthoMatrix="$IsolateAbrv"_orthogroups.tab



echo "$Usage"
echo ""
echo "The following inputs were given:"
echo "MergeHits = $MergeHits"
echo "GoodProts = $GoodProts"
echo "Inflation value = $Inflation"
echo "output will be copied to:"
echo "OutDir = $OutDir"
echo "Files this script will make:"
echo "Config = $Config"
echo "SimilarGenes = $SimilarGenes"
echo "Log_file = $Log_file"
echo "MclInput = $MclInput"
echo "MclOutput = $MclOutput"
echo "OrthoGroups = $OrthoGroups"
echo "OrthoMatrix = $OrthoMatrix"


# ----------------------	Step 2	----------------------
#    Copy the template orthomcl config
#    & Edit fields in this file
#-------------------------------------------------------

# cp /home/armita/testing/armita_orthomcl/orthomcl.config $Config
# cp /home/armita/testing/armita_orthomcl/orthomcl_cluster.config $Config
cp /projects/oldhome/armita/testing/armita_orthomcl/orthomcl_cluster.config $Config

TabName1="$IsolateAbrv"_SimilarSequences
TabName2="$IsolateAbrv"_Ortholog
TabName3="$IsolateAbrv"_InParalog
TabName4="$IsolateAbrv"_CoOrtholog
TabName5="$IsolateAbrv"_interTaxonMatch
sed -i "s/similarSequencesTable=.*/similarSequencesTable="$TabName1"/g" $Config
sed -i "s/orthologTable=.*/orthologTable="$TabName2"/g" $Config
sed -i "s/inParalogTable=.*/inParalogTable="$TabName3"/g" $Config
sed -i "s/coOrthologTable=.*/coOrthologTable="$TabName4"/g" $Config
sed -i "s/interTaxonMatchView=.*/interTaxonMatchView="$TabName5"/g" $Config
sed -i "s/pmatch_cutoff=.*/pmatch_cutoff=50/g" $Config
sed -i "s/evalueExponentCutoff=.*/evalueExponentCutoff=-30/g" $Config

mysql -u armita_orthomcl -parmita_orthomcl -h $HostIP armita_orthomcl -e "drop table if exists BestInterTaxonScore;"
mysql -u armita_orthomcl -parmita_orthomcl -h $HostIP armita_orthomcl -e "drop table if exists BestQueryTaxonScore;"
mysql -u armita_orthomcl -parmita_orthomcl -h $HostIP armita_orthomcl -e "drop table if exists BetterHit;"
mysql -u armita_orthomcl -parmita_orthomcl -h $HostIP armita_orthomcl -e "drop table if exists UniqSimSeqsQueryId;"
mysql -u armita_orthomcl -parmita_orthomcl -h $HostIP armita_orthomcl -e "drop table if exists $TabName1;"
mysql -u armita_orthomcl -parmita_orthomcl -h $HostIP armita_orthomcl -e "drop table if exists $TabName2;"
mysql -u armita_orthomcl -parmita_orthomcl -h $HostIP armita_orthomcl -e "drop table if exists $TabName3;"
mysql -u armita_orthomcl -parmita_orthomcl -h $HostIP armita_orthomcl -e "drop table if exists $TabName4;"
mysql -u armita_orthomcl -parmita_orthomcl -h $HostIP armita_orthomcl -e "drop table if exists $TabName5;"

orthomclInstallSchema $Config install_schema.log


# ----------------------	Step 3	----------------------
#    Run orthoMCL
#         a) parse blast hits
#         b) Load blast results into a database
#         c) Identify pairs of homologous genes
#         d) Write output from the database
#         e) Cluster pairs of homologs into orthogroups
#         f) Parse the OrthoMCL orthogroup output into
#             a matrix that cam be opened in R.
#-------------------------------------------------------

#-- a --
mkdir -p goodProtDir
cp $CurPath/$GoodProts goodProtDir/.
orthomclBlastParser $CurPath/$MergeHits goodProtDir >> $SimilarGenes
#-- b --
ls -lh $SimilarGenes # The database will be 5x the size of this file = ~2.5Gb
orthomclLoadBlast $Config $SimilarGenes
#-- c --
orthomclPairs $Config $Log_file cleanup=yes #<startAfter=TAG>
#-- d --
orthomclDumpPairsFiles $Config
mv mclInput $MclInput
#-- e --
mcl $MclInput --abc -I $Inflation -o $MclOutput
cat $MclOutput | orthomclMclToGroups orthogroup 1 > $OrthoGroups
#-- f --
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
$ProgDir/orthoMCLgroups2tab.py $CurPath/$GoodProts $OrthoGroups > $OrthoMatrix


# ----------------------	Step 4	----------------------
#    Copy orthomcl output into an output directory
#    & cleanup
#-------------------------------------------------------

rm -r goodProtDir
mkdir -p $OutDir
mv $WorkDir/* $OutDir/.
