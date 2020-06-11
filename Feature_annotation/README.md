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
 * D) From ORF fragments - Signal peptide & RxLR motif
 * E) From ORF fragments - Hmm evidence of WY domains
 * F) From ORF fragments - Hmm evidence of RxLR effectors


### A) From Augustus gene models - Identifying secreted proteins

 Required programs:
  * SignalP-4.1
  * TMHMM

 Proteins that were predicted to contain signal peptides were identified using
 the following commands:

screen -a

 ```bash
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
CurPath=$PWD
for Proteome in $(ls gene_pred/N.ditissima/R0905_test/final/final_genes_combined.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/final_genes_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_final_preds
$ProgDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_final_preds_*); do
sbatch $ProgDir/pred_signalP.sh $File signalp-5.0
done
done
```
```bash
for Proteome in $(ls gene_pred/N.ditissima/R0905_test/final/final_genes_combined.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/sub_signalP.sh $Proteome
done

for File in $(ls $SplitDir/*_final_preds_*); do
sbatch $ProgDir/pred_signalP.sh $File signalp-5.0
done
done
```






 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:
 ```bash
 #for Strain in Ag02 Ag05 ND8 R37-15; do
 for Strain in Ag08 P112 Ag06 Ag09_A Ag11_A Ag11_C BGV344 ND9 OPC304 R39-15 R42-15 R68-17 Ag11_B R41-15 R6-17-2 R6-17-3 Ag02 Ag05 ND8; do
 #for Strain in RS305p RS324p; do
	for SplitDir in $(ls -d gene_pred/final_genes_split/N.*/$Strain); do
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

 Some proteins that are incorporated into the cell membrane require secretion.
 Therefore proteins with a transmembrane domain are not likely to represent
 cytoplasmic or apoplastic effectors.

 Proteins containing a transmembrane domain were identified:

 ```bash
 #for Strain in Ag02 Ag05 ND8 R37-15; do
	for Strain in RS305p RS324p; do
 	for Proteome in $(ls gene_pred/codingquary/N.*/$Strain/*/final_genes_combined.pep.fasta); do
 		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
 		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
 		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
 		qsub $ProgDir/submit_TMHMM.sh $Proteome
 	done
done
 ```

 Those proteins with transmembrane domains were removed from lists of Signal peptide containing proteins

 ```bash
 for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt); do
   Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
   Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
   echo "$Organism - $Strain"
   TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
   cat $File | cut -f1 > $TmHeaders
   SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
   OutDir=$(dirname $SigP)
   ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
   cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
 done
 ```
 ```
	N.ditissima - Ag02
	994
	N.ditissima - Ag04
	986
	N.ditissima - Ag05
	996
	N.ditissima - Ag06
	997
	N.ditissima - Ag08
	1058
	N.ditissima - Ag09_A
	1013
	N.ditissima - Ag11_A
	1050
	N.ditissima - Ag11_B
	1019
	N.ditissima - Ag11_C
	1096
	N.ditissima - BGV344
	1000
	N.ditissima - Hg199
	991
	N.ditissima - Hg199_minion
	1033
	N.ditissima - ND8
	996
	N.ditissima - ND9
	1005
	N.ditissima - OPC304
	1016
	N.ditissima - P112
	1029
	N.ditissima - R0905_canu_2017_v2
	965
	N.ditissima - R37-15
	940
	N.ditissima - R39-15
	1019
	N.ditissima - R41-15
	1001
	N.ditissima - R42-15
	1009
	N.ditissima - R45-15
	988
	N.ditissima - R6-17-2
	1005
	N.ditissima - R6-17-3
	989
	N.ditissima - R68-17
	1030
	N.ditissima - RS305p
	1024
	N.ditissima - RS324p
	1007
```

### B) From Augustus gene models - Effector identification using EffectorP

Required programs:

conda install -c bioconda emboss

/usr/lib/x86_64-linux-gnu/libpq.so.5

export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/libpq.so.5
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/libpq.so.5:$LD_LIBRARY_PATH


 * EffectorP.py

```bash
for Proteome in $(ls gene_pred/N.ditissima/R0905_test/final/final_genes_combined.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
BaseName="$Organism"_"$Strain"_EffectorP
OutDir=analysis/effectorP/$Organism/$Strain
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
sbatch $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
done

```

```bash
for Strain in Ag02 Ag05 ND8 R37-15; do
  for File in $(ls analysis/effectorP/*/$Strain/*_EffectorP.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
    cat $File | grep 'Effector' | cut -f1 > $Headers
    Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
    OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
    OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
    cat $OutFileHeaders | wc -l
    Gff=$(ls gene_pred/codingquary/$Organism/$Strain/*/final_genes_appended.gff3)
    EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
    cat $EffectorP_Gff | grep -w 'gene' | wc -l
  done > tmp.txt
done
```

### C) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
#for Strain in Ag02 Ag05 ND8 R37-15; do
	for Strain in RS305p RS324p; do
  for Proteome in $(ls gene_pred/codingquary/N.*/$Strain/*/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/CAZY/$Organism/$Strain
    mkdir -p $OutDir
    Prefix="$Strain"_CAZY
    CazyHmm=dbCAN/dbCAN-fam-HMMs.txt
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
    qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
  done
done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

  ```bash
	#for Strain in Ag02 Ag05 ND8 R37-15; do
	for Strain in RS305p RS324p; do
  for File in $(ls gene_pred/CAZY/N.*/$Strain/*CAZY.out.dm); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $File)
  echo "$Organism - $Strain"
  ProgDir=/home/groups/harrisonlab/dbCAN
  $ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
  CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
  cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
  echo "number of CAZY genes identified:"
  cat $CazyHeaders | wc -l
  Gff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
  CazyGff=$OutDir/"$Strain"_CAZY.gff
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

  SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
  SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
  cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
  CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
  $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
  echo "number of Secreted CAZY genes identified:"
  cat $CazyGffSecreted | grep -w 'gene' | cut -f9 | tr -d 'ID=' | wc -l
  done
done
```
  N.ditissima - Ag04
  number of CAZY genes identified:
  726
  number of Secreted CAZY genes identified:
  283
  N.ditissima - R45-15
  number of CAZY genes identified:
  722
  number of Secreted CAZY genes identified:
  287
  N.ditissima - Hg199
  number of CAZY genes identified:
  721
  number of Secreted CAZY genes identified:
  285
  N.ditissima - R0905_canu_2017_v2
  number of CAZY genes identified:
  719
  number of Secreted CAZY genes identified:
  286

	N.ditissima - Ag02
	number of CAZY genes identified:
	720
	number of Secreted CAZY genes identified:
	281
	N.ditissima - Ag05
	number of CAZY genes identified:
	722
	number of Secreted CAZY genes identified:
	279
	N.ditissima - ND8
	number of CAZY genes identified:
	717
	number of Secreted CAZY genes identified:
	275
	N.ditissima - R37-15
	number of CAZY genes identified:
	713
	number of Secreted CAZY genes identified:
	272
	N.ditissima - RS305p
	number of CAZY genes identified:
	721
	number of Secreted CAZY genes identified:
	287
	N.ditissima - RS324p
	number of CAZY genes identified:
	721
	number of Secreted CAZY genes identified:
	279

Note - the CAZY genes identified may need further filtering based on e value and
cuttoff length - see below:

Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage

* For fungi, use E-value < 1e-17 and coverage > 0.45

* The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)

### Alignment of raw reads vs the Fus2 genome

Sequence data for isolates with a data from a single sequencing run was aligned against the R0905 genome

```bash
  Reference=$(ls repeat_masked/N.ditissima/Ref_Genomes/R0905_canu_2017_v2/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  for StrainPath in $(ls -d qc_dna/paired/N.ditissima/new_isolates_2017/*); do
    Organism=$(echo $StrainPath | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_R0905_canu_2017/
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
  done
  ```

  ## D) Secondary metabolites (Antismash and SMURF)

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org


Results of web-annotation of gene clusters within the assembly were downloaded to
the following directories:

```bash
  for Assembly in $(ls repeat_masked/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
    mkdir -p $OutDir
  done
```

```bash
  for Zip in $(ls analysis/secondary_metabolites/antismash/*/*/*.zip); do
    OutDir=$(dirname $Zip)
    unzip -d $OutDir $Zip
  done
```

```bash
  for AntiSmash in $(ls analysis/secondary_metabolites/antismash/*/*/*/*.final.gbk); do
    Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
    Prefix=$OutDir/WT_antismash
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/antismash2gff.py --inp_antismash $AntiSmash --out_prefix $Prefix

    # Identify secondary metabolites within predicted clusters
    printf "Number of secondary metabolite detected:\t"
    cat "$Prefix"_secmet_clusters.gff | wc -l
    GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
    bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
    cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_antismash_secmet_genes.txt
    bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
    printf "Number of predicted proteins in secondary metabolite clusters:\t"
    cat "$Prefix"_secmet_genes.txt | wc -l
    printf "Number of predicted genes in secondary metabolite clusters:\t"
    cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l

      # Identify cluster finder additional non-secondary metabolite clusters
      printf "Number of cluster finder non-SecMet clusters detected:\t"
      cat "$Prefix"_clusterfinder_clusters.gff | wc -l
      GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
      bedtools intersect -u -a $GeneGff -b "$Prefix"_clusterfinder_clusters.gff > "$Prefix"_clusterfinder_genes.gff
      cat "$Prefix"_clusterfinder_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_clusterfinder_genes.txt

      printf "Number of predicted proteins in cluster finder non-SecMet clusters:\t"
      cat "$Prefix"_clusterfinder_genes.txt | wc -l
      printf "Number of predicted genes in cluster finder non-SecMet clusters:\t"
      cat "$Prefix"_clusterfinder_genes.gff | grep -w 'gene' | wc -l
  done
```

These clusters represented the following genes. Note that these numbers just
show the number of intersected genes with gff clusters and are not confirmed by
function

```
F.venenatum - WT
Number of secondary metabolite detected:	35
Number of predicted proteins in secondary metabolite clusters:	986
Number of predicted genes in secondary metabolite clusters:	977
Number of cluster finder non-SecMet clusters detected:	86
Number of predicted proteins in cluster finder non-SecMet clusters:	2829
Number of predicted genes in cluster finder non-SecMet clusters:	2813
```

SMURF was also run to identify secondary metabolite gene clusters.

Genes needed to be parsed into a specific tsv format prior to submission on the
SMURF webserver.

```bash
  Gff=$(ls gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3)
  OutDir=analysis/secondary_metabolites/smurf/F.venenatum/WT
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/gff2smurf.py --gff $Gff > $OutDir/WT_genes_smurf.tsv
```

SMURF output was received by email and downloaded to the cluster in the output
directory above.

Output files were parsed into gff format:

```bash
  OutDir=analysis/secondary_metabolites/smurf/F.venenatum/WT
  Prefix="WT"
  GeneGff=gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3
  SmurfClusters=$OutDir/Secondary-Metabolite-Clusters.txt
  SmurfBackbone=$OutDir/Backbone-genes.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/smurf2gff.py --smurf_clusters $SmurfClusters --smurf_backbone $SmurfBackbone > $OutDir/Smurf_clusters.gff
  bedtools intersect -wo -a $GeneGff -b $OutDir/Smurf_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > $OutDir/"$Prefix"_smurf_secmet_genes.tsv
```

Total number of secondary metabolite clusters:

```bash
for Assembly in $(ls repeat_masked/*/*/illumina_assembly_ncbi/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
mkdir -p $OutDir
GeneGff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
AntismashClusters=$(ls analysis/secondary_metabolites/antismash/$Organism/$Strain/*_secmet_clusters.gff)
SmurfClusters=$(ls analysis/secondary_metabolites/smurf/$Organism/$Strain/Smurf_clusters.gff)
echo "Total number of Antismash clusters"
cat $AntismashClusters | wc -l
echo "Total number of SMURF clusters"
cat $SmurfClusters | wc -l
echo "number of Antismash clusters intersecting Smurf clusters"
bedtools intersect -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of Antismash clusters not intersecting Smurf clusters"
bedtools intersect -v -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of smurf clusters intersecting antismash clusters"
bedtools intersect -a $SmurfClusters -b $AntismashClusters | wc -l
echo "number of smurf clusters not intersecting antismash clusters"
bedtools intersect -v -a $SmurfClusters -b $AntismashClusters | wc -l
done
```