# Cassis 

### Requirements

```bash
conda activate meme_v4
# Cassis is only compatible with version 4 of MEME.
conda install meme=4.12.0
# Note: For version 5 of MEME, the output of FIMO will be called fimo.tsv and this need to be modified in cassis.pl

# Cassus is intalled in /scratch. Add the following line to your profile
PATH=${PATH}:/scratch/software/cassis/CASSIS_Linux_64bit
```

### Typical run

```bash
screen -a

ProjDir=$(ls -d /path/to/directory)
WorkDir=$ProjDir/analysis/promoters/cassis/all_genes # If all secmet genes are used
mkdir -p $WorkDir
cd $WorkDir

AnnotTab=$(ls $ProjDir/analysis/annotation_tables/$Organism/$Strain/*_gene_table.tsv) # Annotation tables 
Assembly=$(ls $ProjDir/repeat_masked/$Organism/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa) # Genome assembly    
Genes=$(ls $ProjDir/gene_pred/codingquarry/$Organism/$Strain/final/final_genes_appended_renamed.gff3) # Final gff3 file
Interpro=$(ls $ProjDir/gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv) # Interproscan annotation tsv

# Extract gene name, contigs and coordinates
cat $Genes | grep 'mRNA' | sed 's/ID=//g' | sed "s/;.*//g" | awk '{ print $9 "\t" $1 "\t" $4 "\t" $5 "\t" $7}' > cassis.tsv
CassisTSV=cassis.tsv

# Run cassis

for Cluster in $(cat $AnnotTab | cut -f13 | grep 'contig' | sort -n -k3 -t'_' | sed 's/;.*//p' | uniq); do
    echo $Cluster
    mkdir $WorkDir/$Cluster
    cat $AnnotTab | cut -f1,13 | grep -w "$Cluster" | cut -f1 | grep '.t1' > $WorkDir/$Cluster/headers.txt
for GeneID in $(cat $WorkDir/$Cluster/headers.txt); do
    echo $GeneID
    mkdir -p $WorkDir/$Cluster/$GeneID
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Promoter_analysis
    OutDir=$WorkDir/$Cluster
    Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
        while [ $Jobs -gt 60 ]; do
            sleep 5m
            printf "."
            Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
        done
    printf "\n"
    sbatch $ProgDir/cassis.sh $Assembly $CassisTSV $GeneID $OutDir
done
done
```






  

```bash
ProjDir=$(ls -d /projects/fusarium_venenatum)
cd $ProjDir
for Cluster in $(ls -d analysis/promoters/cassis/all_genes/contig* | rev | cut -f1 -d '/' | rev | sort -n -k3 -t'_'); do
ClusterDir=$(ls -d analysis/promoters/cassis/all_genes/${Cluster})
echo ""
for Results in $(ls $ClusterDir/*/*_log.txt); do
Anchor=$(echo $Results | rev | cut -f2 -d '/' | rev)
# if [[ $Cluster == "" ]]; then
# Cluster="NA"
# fi
if $(grep -q 'No cluster prediction' $Results); then
printf "${Cluster}\t${Anchor}\tNA\tNA\n"
elif grep 'Computing CLUSTER PREDICTIONS' $Results; then
Best=$(cat $Results | grep -A2 '(7) Computing CLUSTER PREDICTIONS' | tail -n1 | sed -r "s&^\s+&&g" | cut -f1 -d ' ')
Fimo=$(ls $ClusterDir/$Anchor/$Anchor/fimo/$Best/fimo.txt)
Motif=$(cat $Fimo | head -n2 | tail -n1 | cut -f1)
printf "${Cluster}\t${Anchor}\t${Best}\t${Motif}\n"
else
printf "${Cluster}\t${Anchor}\tNA\tNA\n"
fi
done | grep -v 'CLUSTER PREDICTIONS' | grep -v ':('
done > analysis/promoters/cassis/all_genes/cassis_summary.tsv
```