# Gene feature annotation tools

1. Interproscan. Provides functional analysis of proteins by classifying them into families and predicting domains.

2. Swissprot

3. SignalP

4. EffectorP

7. Antismash

## 1.Interproscan

Interproscan was used to give gene models functional annotations.

### Requirements

```bash
# The lastest versions of interproscan are installed in /data/scratch/gomeza/prog/Interproscan and used by in the run_interproscan.sh script.
# The current version of interproscan only works with Java version 11. The next line has to be added to your profile.
PATH=/data/scratch/gomeza/prog/java/jdk-11.0.4/bin:${PATH}

. ~/.profile # Refresh your profile
```

```bash
# This command will split your gene fasta file and run multiple interproscan jobs.
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  for Genes in $(ls path/to/the/final_genes_appended_renamed.pep.fasta); do
    echo $Genes
    $ProgDir/interproscan.sh $Genes
  done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following commands:

```bash
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  for Proteins in $(ls path/to/the/final_genes_appended_renamed.pep.fasta); do
    Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    echo $Strain
    InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
    $ProgDir/append_interpro.sh $Proteins $InterProRaw
  done
```


## 2.SwissProt


### Requirements

```bash
# No requirements to run Swissprot
# Uniprot databases are downloaded to /projects/dbUniprot

# Intructions to create a database if needed.
dbFasta=$(ls /projects/dbUniprot/swissprot_2020_June/uniprot_sprot.fasta)
dbType="prot"
Prefix="uniprot_sprot"
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix 
```

### Typical run


```bash
for Proteome in $(ls gene_pred/N.ditissima/R0905_test/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../dbUniprot/swissprot_2020_June
SwissDbName=uniprot_sprot
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```


## Effector genes

Putative pathogenicity and effector related genes were identified within Braker
gene models using a number of approaches:

 * A) From Augustus gene models - Identifying secreted proteins
 * B) From Augustus gene models - Effector identification using EffectorP


## 3. SignalP and TMHMM

From Augustus gene models - Identifying secreted proteins

### Requirements

```
# SignalP and TMHMM needed. Add these path to your profile (if not done before)

PATH=${PATH}:/data/scratch/gomeza/prog/signalp/signalp-5.0b/bin
PATH=${PATH}:/data/scratch/gomeza/prog/signalp/signalp-4.1
PATH=${PATH}:/data/scratch/gomeza/prog/tmhmm-2.0c/bin

. ~/.profile # Refresh your profile
```

Proteins that were predicted to contain signal peptides were identified using signalP

```bash
for Strain in Strain1 Strain2 Strain3; do # Add your strains name
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  CurPath=$PWD
  for Proteome in $(ls gene_pred/*/$Strain/final/final_genes_combined.pep.fasta); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  SplitDir=gene_pred/final_genes_split/$Organism/$Strain
  mkdir -p $SplitDir
  BaseName="$Organism""_$Strain"_final_preds
  $ProgDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName # Split your input fasta in 500 genes files
    for File in $(ls $SplitDir/*_final_preds_*); do
    #sbatch $ProgDir/pred_signalP.sh $File signalp 
    #sbatch $ProgDir/pred_signalP.sh $File signalp-3.0 # Recommended for oomycetes
    sbatch $ProgDir/pred_signalP.sh $File signalp-4.1 # Recommended for fungi
    #sbatch $ProgDir/pred_signalP.sh $File signalp-5.0
    done
  done
done
```

 The batch files of predicted secreted proteins needed to be combined into a single file for each strain. This was done with the following commands:

 ```bash
 for Strain in Strain1 Strain2 Strain3; do  # Add your strains name
	for SplitDir in $(ls -d gene_pred/final_genes_split/*/$Strain); do
		Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
		Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
		InStringAA=''
		InStringNeg=''
		InStringTab=''
		InStringTxt=''
		SigpDir=final_genes_signalp-4.1
		for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do
			InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";
			InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";
			InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
			InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";
		done
		cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
		cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
		tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
		cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
	done
done
 ```
Some proteins that are incorporated into the cell membrane require secretion. Therefore proteins with a transmembrane domain are not likely to represent cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

 ```bash
for Strain in Strain1 Strain2 Strain3; do  # Add your strains name
 	for Proteome in $(ls gene_pred/codingquary/*/$Strain/final/final_genes_combined.pep.fasta); do
 		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
 		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
 		ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
 		qsub $ProgDir/TMHMM.sh $Proteome
 	done
done
 ```

 Those proteins with transmembrane domains were removed from lists of Signal peptide containing proteins

 ```bash
for File in $(ls gene_pred/trans_mem/$Organism/$Strain/*_TM_genes_neg.txt); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
  cat $File | cut -f1 > $TmHeaders
  SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
  OutDir=$(dirname $SigP)
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  $ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
  cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
done
```


## 4. EffectorP

From Augustus gene models - Effector identification using EffectorP

### Requirements
```
# This line need to added to profile or used to execute EffectorP version 2.0 directly
PATH=${PATH}:/scratch/software/EffectorP-2.0/Scripts
```
```bash
for Proteome in $(ls path/to/final/final_genes_appended_renamed.pep.fasta); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  BaseName="$Organism"_"$Strain"_EffectorP
  OutDir=analysis/effectorP/$Organism/$Strain
  Version=v2 # Version 2.0 or 3.0
  OutDir=analysis/effectorP/$Version/$Organism/$Strain
  #EffectorP.py -o "$BaseName".txt -E "$BaseName".fa -i $Proteome
  #mv "$BaseName".txt $OutDir
  #mv "$BaseName".fa $OutDir
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  sbatch $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir $Version
done
```


```bash
for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
  cat $File | grep 'Effector' | cut -f1 > $Headers
  Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
  OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
  OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
  cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
  cat $OutFileHeaders | wc -l
  Gff=$(ls gene_pred/codingquarry_cuff_final/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3)
  EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
  $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
  cat $EffectorP_Gff | grep -w 'gene' | wc -l
done > tmp.txt
```

## 5. Identification of MIMP-flanking genes

Miniature impala (mimp) sequeces are found in promotor regions of SIX genes in fusarium.

```bash
  for Assembly in $(ls path/to/*_contigs_unmasked.fa); do
    Organism=$(echo "$Assembly" | rev | cut -d '/' -f4 | rev)
    Strain=$(echo "$Assembly" | rev | cut -d '/' -f3 | rev)
    GeneGff=$(ls gene_pred/final_genes/$Organism/"$Strain"/final/final_genes_appended_renamed.gff3)
    OutDir=analysis/mimps/$Organism/$Strain
    mkdir -p "$OutDir"
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
    $ProgDir/mimp_finder.pl $Assembly $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
    $ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
    echo "The number of mimps identified:"
    cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
    bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
    echo "The following transcripts intersect mimps:"
    MimpProtsTxt=$OutDir/"$Strain"_prots_in_2kb_mimp.txt
    MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
    cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq > $MimpProtsTxt
    cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | cut -f1 -d '.'| sort | uniq > $MimpGenesTxt
    cat $MimpProtsTxt | wc -l
    cat $MimpGenesTxt | wc -l
    echo ""
  done
```


## 6. CAZY proteins

Carbohydrte active enzymes were identified from the CAZy database

```bash
for Strain in Strain1 Strain2; do # List of isolates
  for Proteome in $(ls path/to/pep/fasta/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/CAZY/$Organism/$Strain
    mkdir -p $OutDir
    Prefix="$Strain"_CAZY
    CazyHmm=../dbCAN/dbCAN-HMMdb-V8.txt # databases are in /projects
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
    sbatch $ProgDir/hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
  done
done
```
```bash
 for File in $(ls path/to/*CAZY.out.dm); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $File)
  echo "$Organism - $Strain"
  ProgDir=/projects/dbCAN # Script from dbCAN tools
  $ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps # Creates a file with CAZy module and gene

  CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
  cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders # Extract gene names
  echo "Number of CAZY genes identified:"
  cat $CazyHeaders | wc -l

  Gff=$(ls path/to/final/gff3/file/final_genes_appended_renamed.gff3)
  CazyGff=$OutDir/"$Strain"_CAZY.gff
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  $ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff # Creates a gff for all CAZymes
  
  SecretedProts=$(ls path/to/secreted/proteins/*signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
  SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
  cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
  CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
  $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted # Creates a gff for secreted CAZymes
  echo "Number of Secreted CAZY genes identified:"
  cat $CazyGffSecreted | grep -w 'gene' | cut -f9 | tr -d 'ID=' | wc -l
  done
  ```

## 7. Antismash

Antismash was run to identify clusters of secondary metabolite genes within the genome. Antismash was run using the webserver at: http://antismash.secondarymetabolites.org.

Results of web-annotation of gene clusters within the assembly were downloaded to the following directories

```bash
for AntiSmash in $(ls path/to/antismash/output/gbk/file/*appended.gbk); do
  Organism=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $AntiSmash | rev | cut -f2 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
  Prefix=$OutDir/"Strain"_antismash_results
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  # Antismash v5 output to gff file
  $ProgDir/antismash2gffv5.py --inp_antismash $AntiSmash --out_prefix $Prefix 
  #$ProgDir/antismash2gff.py --inp_antismash $AntiSmash --out_prefix $Prefix # Use only for antismash v4.2 output
  printf "Number of secondary metabolite detected:\t"
  cat "$Prefix"_secmet_clusters.gff | wc -l
  GeneGff=path/to/final/gff3/file/final_genes_appended_renamed.gff3
  bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
  cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_antismash_secmet_genes.txt
  bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -e "s/;Parent=g\w+//g" | perl -p -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
  printf "Number of predicted proteins in secondary metabolite clusters:\t"
  cat "$Prefix"_secmet_genes.txt | wc -l
  printf "Number of predicted genes in secondary metabolite clusters:\t"
  cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l
done

# Antismash output correction. Some gene names contain ;. Remove manually with the following command.
# First sed command removes ;. Second and Third remove the cluster kind information (optional)
cat "Strain"_antismash_results_secmet_genes.tsv | sed 's/;//p' | sed 's/;.*//p' | sed 's/Kin.*//p' > "Strain"_antismash_results_secmet_genes_corrected.tsv
 ```

## 8. Looking for Transcription factors

```bash
  for Interpro in $(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv); do
    Organism=$(echo $Interpro | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Interpro | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/transcription_factors/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
    $ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
    echo "total number of transcription factors"
    cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
    cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l
    # Gene ID rather than transcript ID
    cat $OutDir/"$Strain"_TF_gene_headers.txt | sed -e "s/.t.*//g" > $OutDir/"$Strain"_TF_geneid_headers.txt
  done
```
