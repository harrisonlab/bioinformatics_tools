### Align mycelium reads to Hg199 assembly

```bash
for Assembly in $(ls assembly/miniasm/N.ditissima/Hg199/racon_10/pilon/filtered_contigs/repeat_masked/Hg199_miniasm_contigs_unmasked.fa)
do
Strain=Hg199
Organism=N.ditissima
echo "$Organism - $Strain"
for FileF in $(ls qc_rna/RNAseq/N.ditissima/Hg199/mycelium/F/*_trim.fq.gz)
do
FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_trim/_2_trim/g')
echo $FileF
echo $FileR
Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
echo "$Timepoint"
Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
OutDir=alignment/star/$Organism/$Strain/vs_miniasm/$Timepoint/$Sample_Name
Preindexlength=12 # Depending genome size. E.g 45Mb = 12
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir $Preindexlength
done
done
```

### Align all timepoints to apple genome

```bash
for Assembly in $(ls apple_genome/GDDH13_1-1_formatted.fasta)
do
Strain=GD
Organism=M.domestica
echo "$Organism - $Strain"
for FileF in $(ls qc_rna/RNAseq/N.ditissima/M9/t0/F/*_trim.fq.gz)
do
FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_trim/_2_trim/g')
echo $FileF
echo $FileR
Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
echo "$Timepoint"
Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
Preindexlength=13 # Depending genome size. E.g 45Mb = 12, Apple=13
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir $Preindexlength
done
done
```

### Gzip output files to save space on the disk and allow star to run correctly downstream. ONLY RUN THIS ONCE

```bash
for AlignDir in $(ls -d alignment/star/N.ditissima/*/*/*)
do
    cat $AlignDir/star_aligmentUnmapped.out.mate1 | gzip -cf > $AlignDir/star_aligmentUnmapped.out.mate1.fq.gz
    cat $AlignDir/star_aligmentUnmapped.out.mate2 | gzip -cf > $AlignDir/star_aligmentUnmapped.out.mate2.fq.gz
done
```

### Align unmapped reads to Hg199 genome

```bash
for Assembly in $(ls apple_genome/GDDH13_1-1_formatted.fasta)
do
Strain=GD
Organism=M.domestica
echo "$Organism - $Strain"
for FileF in $(ls qc_rna/RNAseq/N.ditissima/M9/t0/F/*_trim.fq.gz)
do
FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_trim/_2_trim/g')
echo $FileF
echo $FileR
Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
echo "$Timepoint"
Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
Preindexlength=13 # Depending genome size. E.g 45Mb = 12, Apple=13
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir $Preindexlength
done
done
```

# Feature counts

```bash
for BamFile in $(ls -d alignment/star/M.domestica/GD/t1/GD_4A4/star_aligmentAligned.sortedByCoord.out.bam); do
Gff=apple_genome/gene_models_20170612.gff3
Time=$(echo $BamFile | rev | cut -f3 -d '/' | rev)
Prefix=$(echo $BamFile | rev | cut -f2 -d '/' | rev)
OutDir=alignment/star/M.domestica/GD/$Time/$Prefix/featureCounts
echo $Prefix
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
sbatch $ProgDir/featurecounts.sh $BamFile $Gff $OutDir $Prefix
done
```


cut -f1,7 counts018347.txt > tmp1
cut -f7 counts019035.txt > tmp2
cut -f7 SRR074122.txt > tmp3
cut ...
paste tmp1 tmp2 tmp3 tmp... > summary.txt


```bash
mkdir -p analysis/transposons
cat gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv | grep -e 'IPR000477' -e 'IPR012337' -e 'IPR018289' -e 'PF03221' -e 'PF00078' -e 'IPR025476' -e 'IPR008906' -e 'transpos' -e 'integrase' | cut -f1 | sort | uniq > analysis/transposons/transposon_headers.txt
```

```bash
for Alignment in $(ls alignment/star/M.domestica/GD/t1/GD_4A4/star_aligmentAligned.sortedByCoord.out.bam); do
Gff=apple_genome/gene_models_20170612.gff3
OutDir=$(dirname $Alignment)
mkdir -p $OutDir/fpkm
cufflinks -p 8 -o $OutDir/fpkm -G $Gff $Alignment
done

fpkm_files=$(ls norm/v2/genes.fpkm_tracking | sed -r 's/\n/ /g')
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
#Cazy=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt
#EffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_
Cazy=/projects/oldhome/groups/harrisonlab/project_files/fusarium/gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_headers.txt
EffectorP=/projects/oldhome/groups/harrisonlab/project_files/fusarium/analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_headers.txt
#EffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_headers.txt
Transposons=/projects/oldhome/groups/harrisonlab/project_files/fusarium/analysis/transposons/transposon_headers.txt
Metabolites=/projects/oldhome/groups/harrisonlab/project_files/fusarium/analysis/antismash/F.oxysporum_fsp_cepae/Fus2_canu_new/metabolite_cluster_gene_headers.txt
Annotations=/projects/oldhome/groups/harrisonlab/project_files/fusarium/gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
$ProgDir/fpkm_by_gene.py --fpkm_file $fpkm_files --CAZY_headers $Cazy --effectorP_headers $EffectorP  --transposon_headers $Transposons --metabolite_headers $Metabolites --annotation_table $Annotations > analysis/expression/Fus2_expressed_genes_PDB.tsv
```





```bash
mkdir -p analysis/transposons
cat /data/scratch/gomeza/gene_pred_vAG/interproscan/Ref_Genomes/N.ditissima/Hg199/Hg199_interproscan.tsv | grep -e 'IPR000477' -e 'IPR012337' -e 'IPR018289' -e 'PF03221' -e 'PF00078' -e 'IPR025476' -e 'IPR008906' -e 'transpos' -e 'integrase' | cut -f1 | sort | uniq > analysis/transposons/transposon_headers.txt
```

for Alignment in $(ls alignment/star/N.ditissima/Hg199/vs_miniasm/mycelium/Hg199_1/star_aligmentAligned.sortedByCoord.out.bam); do
Gff=/data/scratch/gomeza/gene_pred_vAG/codingquary/Ref_Genomes/N.ditissima/Hg199/final/final_genes_appended_renamed.gff3 
OutDir=$(dirname $Alignment)
mkdir -p $OutDir/fpkm
cufflinks -p 8 -o $OutDir/fpkm -G $Gff $Alignment
done

fpkm_files=$(ls alignment/star/N.ditissima/Hg199/vs_miniasm/mycelium/Hg199_1/fpkm/genes.fpkm_tracking | sed -r 's/\n/ /g')
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
#Cazy=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt
#EffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_
Cazy=/data/scratch/gomeza/gene_pred_vAG/CAZY/Ref_Genomes/N.ditissima/Hg199/Hg199_CAZY_headers.txt
EffectorP=/projects/neonectria_ditissima/analysis_vAG/effectorP/Ref_Genomes/N.ditissima/Hg199/N.ditissima_Hg199_EffectorP_headers.txt
#EffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_headers.txt
Transposons=analysis/transposons/transposon_headers.txt
Metabolites=/projects/neonectria_ditissima/analysis_vAG/secondary_metabolites/antismash/Ref_Genomes/N.ditissima/Hg199/Antismash_v4.2/Hg199_secmet_genes.txt
Annotations=/data/scratch/gomeza/gene_pred_vAG/annotation/Pathogen_vAG/N.ditissima/Hg199/Hg199_gene_table_incl_exp_GROUP.tsv
$ProgDir/fpkm_by_gene.py --fpkm_file $fpkm_files --CAZY_headers $Cazy --effectorP_headers $EffectorP  --transposon_headers $Transposons --metabolite_headers $Metabolites --annotation_table $Annotations > analysis/expression/Nd_expressed_genes.tsv