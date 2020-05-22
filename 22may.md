```bash
Strain=R0905
Organism=N.ditissima
samtools merge -f concatenated/concatenated.bam \
Hg199_1/star_aligmentAligned.sorted.out.bam \
Hg199_2/star_aligmentAligned.sorted.out.bam \
Hg199_3/star_aligmentAligned.sorted.out.bam

Strain=R0905
Organism=N.ditissima
samtools merge -f concatenated/concatenated_unsorted.bam \
Hg199_1/star_aligmentAligned.out.bam \
Hg199_2/star_aligmentAligned.out.bam \
Hg199_3/star_aligmentAligned.out.bam
```

stringtie concatenated.bam -o short_reads.out.gtf -v

stringtie concatenated.bam -o out_concatenated.gtf -p 8 -A gene_abund.tab

sed -i 1,2d out.gtf

```bash
for Assembly in $(ls R0905_contigs_unmasked.fa); do
#Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev) # Edit to set your ouput directory
#Organism=$(echo $Assembly| rev | cut -d '/' -f6 | rev) # Edit to set your ouput directory
#echo "$Organism - $Strain"
OutDir=gene_pred_last
mkdir -p $OutDir
GTF=../out_concatenated.gtf  # GFT file from stringtie/cufflinks output. See Genome-guided_assemblers scripts
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/codingquarry.sh $Assembly $GTF $OutDir
done
```


Akin

```bash
for Assembly in $(ls ../../fusarium_ex_strawberry/repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa)
do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev) # Edit to set your ouput directory
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev) # Edit to set your ouput directory
echo "$Organism - $Strain"
for FileF in $(ls ../qc_rna/RNAseq/N.ditissima/Hg199/mycelium/F/Hg199_1_1_trim.fq.gz)
do
FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_trim/_2_trim/g')
echo $FileF
echo $FileR
Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
echo "$Timepoint"
Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
OutDir=alignment_akin/star/$Organism/$Strain/$Timepoint/$Sample_Name
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir
done
done

for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa)
  do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for FileF in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/control_72hrs_rep1/F/*_trim.fq.gz)
    do
  FileR=../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/control_72hrs_rep1/R/*_trim.fq.gz
    echo $FileF
    echo $FileR
    Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
    echo "$Timepoint"
    Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
    OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
    sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir
  done
done
```

```bash
  BrakerGff=$(ls gene_pred_vAG/braker/Ref_Genomes/*/R0905/N.ditissima_R0905_braker/augustus.gff3)
	Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev)
	Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
	echo "$Organism - $Strain"
	Assembly=$(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
	CodingQuaryGff=gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/out/PredictedPass.gff3
	PGNGff=gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/out/PGN_predictedPass.gff3
	AddDir=gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/additional
	FinalDir=gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/final
	AddGenesList=$AddDir/additional_genes.txt
	AddGenesGff=$AddDir/additional_genes.gff
	FinalGff=$AddDir/combined_genes.gff
	mkdir -p $AddDir
	mkdir -p $FinalDir

	bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
	bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation
	$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
	$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/codingquary

	$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
	$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
	cp $BrakerGff $FinalDir/final_genes_Braker.gff3
	$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
	cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
	cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
	cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
	cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

	GffBraker=$FinalDir/final_genes_CodingQuary.gff3
	GffQuary=$FinalDir/final_genes_Braker.gff3
	GffAppended=$FinalDir/final_genes_appended.gff3
	cat $GffBraker $GffQuary > $GffAppended
    ```

    ```bash
for Assembly in $(ls ../../../../../../../../assembly/canu/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev) # Edit to set your ouput directory
Organism=$(echo $Assembly | rev | cut -d '/' -f6 | rev) # Edit to set your ouput directory
echo "$Organism - $Strain"
OutDir=gene_pred/braker/$Organism/$Strain
AcceptedHits=../concatenated.bam # Concatenatented alignment files can be used
GeneModelName="$Organism"_"$Strain"_braker 
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```
