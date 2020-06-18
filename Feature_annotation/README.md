# Gene feature annotation tools

1. Interproscan. Provides functional analysis of proteins by classifying them into families and predicting domains.

2. Swissprot

3. SignalP

4. EffectorP

## Interproscan

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






## B) SwissProt

```bash
for Proteome in $(ls gene_pred/N.ditissima/R0905_test/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=/projects/dbUniprot/swissprot_2020_June
SwissDbName=uniprot_sprot
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done

```

```bash
	for SwissTable in $(ls gene_pred/swissprot/*/*/*_hits.tbl); do
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_v2017_tophit_parsed.tbl
		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
	done
```

## Effector genes

Putative pathogenicity and effector related genes were identified within Braker
gene models using a number of approaches:

 * A) From Augustus gene models - Identifying secreted proteins
 * B) From Augustus gene models - Effector identification using EffectorP


### A) From Augustus gene models - Identifying secreted proteins

### Requirements

```
# SignalP and TMHMM needed. Add these path to your profile (if not done before)

PATH=${PATH}:/data/scratch/gomeza/prog/signalp/signalp-5.0b/bin
PATH=${PATH}:/data/scratch/gomeza/prog/signalp/signalp-4.1
PATH=${PATH}:/data/scratch/gomeza/prog/tmhmm-2.0c/bin

. ~/.profile # Refresh your profile
```

Proteins that were predicted to contain signal peptides were identified using the following commands:

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

 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

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



## B) From Augustus gene models - Effector identification using EffectorP

### Requirements
```
# This line need to added to profile or used to execute EffectorP directly
PATH=${PATH}:/scratch/software/EffectorP-2.0/Scripts
```
```bash
for Proteome in $(ls path/to/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
BaseName="$Organism"_"$Strain"_EffectorP
OutDir=analysis/effectorP/$Organism/$Strain
EffectorP.py -o "$BaseName".txt -E "$BaseName".fa -i $Proteome
mv "$BaseName".txt $OutDir
mv "$BaseName".fa $OutDir
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
#sbatch $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
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
