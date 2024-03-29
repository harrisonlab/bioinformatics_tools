# Gene prediction tools

Tools used in the prediction and annotation of genes

1. BUSCO: Indentification of single-copy orthologs that should be highly conserved among the closely related species. This tool can be used to evaluated genome completeness and for phylogenetic analysis.

2. Braker: A pipeline for automated prediction of protein coding genes using GeneMark-ES/ET and AUGUSTUS. This allow ab initio and gene model training predictions.

3. Cufflinks: Assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples.

4. Stringtie: Fast and high efficient assembler of RNA-Seq alignments (e.g STAR output) into transcripts.

5. CodinQuarry: Generalised hidden Markov models gene prediction tool. This tool is highly accurate predicting fungal genes, using "pathogen mode".

Additional tools in this folder

- add_CodingQuarry_features.pl
- gene_list_to_gff.pl
- gff_rename_genes.py
- gff2fasta.pl
- remove_dup_features.py

## Busco

Indentification of Benchmarking Universal Single-Copy Orthologs in the genome usgin BUSCO databases.

### Requirements

```bash
# The latest version of BUSCO might be incompatible with some dependencies already installed in your env.
# Create a new env for BUSCO.
conda create -n BUSCO -c bioconda -c conda-forge busco=4.0.6 # This takes time
conda activate BUSCO

# The latest BUSCO databases, odb10, are in /projects/dbBusco/
# Busco databases odb9 are in /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9
```

### Typical run

```bash
  for Assembly in $(ls path/to/genome/assembly/*.fasta); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev) # Edit to set your ouput directory
    Organism=$(echo $Assembly | rev | cut -d '/' -f6 | rev) # Edit to set your ouput directory
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
    sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
  done
```


## Braker

Braker is a combination of GeneMark-ET and Augustus used for gene annotation.

### Requirements

Braker requires the instalation of multiple perl libraries, Augustus and Genemark. There are conda packages for Augustus and Braker. However, these packages are often lagging behing the releasess of AUGUSTUS and BRAKER. Therefore manual installation might be necessary.

### Conda installation

```bash
Conda activate gene_pred # e.g.

conda install numpy
conda install -c anaconda boost
conda install -c bioconda bamtools
conda install -c bioconda braker
conda install -c bioconda augustus

# GeneMark-ET will run if a VALID key file resides in your home directory.
cp /home/gomeza/prog/genemark/gm_key_64 ~/.gm_key

### Genemark-ET installation (optional)

# In case you need/want to install GeneMark in your home directory. If so, braker_fungi.sh script needs to be edited with the new Genemark path.

# Download GeneMark-ES/ET/EP ver 4.57_lic from http://exon.gatech.edu/GeneMark/license_download.cgi. to your home directory, e.g. /home/gomeza/prog/genemark

tar -xvf gmes_linux_64.tar.gz
gunzip gm_key_64.gz

cd #change to home directory
cp /path/to/genemark/gm_key_64 ~/.gm_key  

### Important configuration

# Current version of AUGUSTUS and Braker can be installed using conda packages (braker version 2.1.5 and AUGUSTUS 3.3.3). However, some perl scripts in the conda package of AUGUSTUS are incompatible with this version of braker, producing unexpected output files. These were updated but not included in the conda package yet. Therefore, replace them in your conda installation before running braker.

cp /home/gomeza/miniconda3/envs/gene_pred/bin/filterGenesIn_mRNAname.pl /home/USER_ID/miniconda3/USER_ENV/bin
```

### Typical run

```bash
  for Assembly in $(ls path/to/softmasked/assembly/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev) # Edit to set your ouput directory
    Organism=$(echo $Assembly | rev | cut -d '/' -f6 | rev) # Edit to set your ouput directory
    echo "$Organism - $Strain"
    OutDir=gene_pred/braker/$Organism/$Strain
    AcceptedHits=path/to/your/spliced/aligments/files/*_aligmentAligned.sortedByCoord.out.bam # STAR output, see Genome_aligners folder
    #AcceptedHits=alignment/concatenated.bam # Concatenatented alignment files can be used
    GeneModelName="$Organism"_"$Strain"_braker 
    ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Gene_prediction
    sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```


## Transcriptome assembly

Cufflinks is older and has not been updated in years. It accepts aligned RNA-Seq reads.

Stringtie is faster but might struggle with large alignment files.

### Requirements

```bash
# Conda installation for stringtie

conda install stringtie

# Cufflinks is installed in the oldhome directory. 

nano ~/.profile # Edit profile
PATH=${PATH}:/projects/oldhome/armita/prog/cufflinks/cufflinks-2.2.1.Linux_x86_64 # Add to profile and save
. ~/.profile # Refresh your profile
```


### Cufflinks typical run

```bash
  for Assembly in $(ls path/to/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) # Edit to set your ouput directory
    Organism=$(echo $Assembly| rev | cut -d '/' -f3 | rev) # Edit to set your ouput directory
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
    mkdir -p $OutDir
    AcceptedHits=path/to/your/spliced/aligments/files.bam
    ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Gene_prediction
    sbatch $ProgDir/cufflinks.sh $AcceptedHits $OutDir
  done
```

### Stringtie typical run

```bash
  for Assembly in $(ls path/to/unmasked/genome/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev) # Edit to set your ouput directory
    Organism=$(echo $Assembly| rev | cut -d '/' -f6 | rev) # Edit to set your ouput directory
    echo "$Organism - $Strain"
    OutDir=gene_pred/stringtie/$Organism/$Strain/concatenated_prelim
    mkdir -p $OutDir
    AcceptedHits=path/to/your/spliced/aligments/files.bam
    ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Gene_prediction
    qbatch $ProgDir/stringties.sh $AcceptedHits $OutDir
   done
```

## CodingQuarry

CodingQuarry in pathogen mode is used to predict aditional genes and added to braker predictions

### Requirements

```bash
# 1st - Conda installation
# CodingQuarry requires a conda environment with python 2.7
# e.g. conda create --name python_py27 python=2.7
mamba install codingquarry

# 2nd - The environmental variable QUARRY_PATH is set in your profile (needed for CodingQuarry)
nano ~/.bash_profile # Edit profile
export QUARRY_PATH="/mnt/shared/scratch/"USER_ID"/apps/conda/envs/"USER_ENV_py27"/opt/codingquarry-2.0/QuarryFiles" # Add to profile and save

# 3rd - Edit python shebangs of CufflinksGTF_to_CodingQuarryGFF3.py and Quarry scripts to run on py27 (e.g. /mnt/shared/scratch/agomez/apps/conda/envs/python_py27/bin/python2.7)

# 4th - run_CQ-PM_stranded.sh and run_CQ-PM_unstranded.sh scripts included in codingquarry are used to run CQ pathogen mode using signalp 4.1. 
# The script in this folder was edited to use signalp5 installed on /home/agomez/scratch/apps/prog/signalp-5.0b
```

### Typical run

Note: run_CQ-PM_stranded.sh and run_CQ-PM_unstranded.sh scripts are included in cndigquarry scripts are used to run CQ pathogen mode using signalp 4.1. The script in this folder was edited to use signalp5. 

```bash
  for Assembly in $(ls path/to/unmasked/genome/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev) # Edit to set your ouput directory
    Organism=$(echo $Assembly| rev | cut -d '/' -f6 | rev) # Edit to set your ouput directory
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain/
    mkdir -p $OutDir
    GTF=path/to/RNAseq/alignment/assembly/*.gtf # GFT file from stringtie/cufflinks output. See Genome-guided_assemblers scripts
    ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Gene_prediction
    sbatch $ProgDir/codingquarry.sh $Assembly $GTF $OutDir
  done
```

### Add additional transcripts to Braker gene models.


Additional transcripts predicted by CodingQuarry are added to the final gene models.

```bash
  # The following perl scripts requires the installation of some libraries. Run these commands in a perly environment.
  # Install the required libraries (if any) using cpanm
  # mamba install perl-app-cpanminus
  # mamba install perl-bioperl
  # cpanm --force Bio::DB::GenBank 
  # cpanm --force Bio::Perl

  BrakerGff=$(ls path/to/braker/gene/models/augustus.hints.gff3)
	Strain=$(echo $BrakerGff| rev | cut -d '/' -f2 | rev)
	Organism=$(echo $BrakerGff | rev | cut -d '/' -f3 | rev)
	echo "$Organism - $Strain"
	Assembly=$(ls path/to/softmasked/genome/assembly/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
	CodingQuarryGff=path/to/codingquarry/gff3/$Organism/$Strain/out/PredictedPass.gff3
	PGNGff=path/to/codingquarry/pathogen/mode/gff3/$Organism/$Strain/out/PGN_predictedPass.gff3
	AddDir=gene_pred/codingquary/$Organism/$Strain/additional # Additional transcripts directory
	FinalDir=gene_pred/codingquary/$Organism/$Strain/final # Final directory
	AddGenesList=$AddDir/additional_genes.txt
	AddGenesGff=$AddDir/additional_genes.gff
	FinalGff=$AddDir/combined_genes.gff
	mkdir -p $AddDir
	mkdir -p $FinalDir

  # Create a list with the additional transcripts in CondingQuarry gff (and CQPM) vs Braker gene models
	bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
	bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
  
  # Creat Gff file with the additional transcripts
	ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Gene_prediction
	$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
	$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
	
  # Create a final Gff file with gene features
	$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3

  # Create fasta files from each gene feature in the CodingQuarry gff3
	$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary

  # Create fasta files from each gene feature in the Braker gff3
  # If any "Possible precedence issue with control flow operator" error, change "or" for "||" operator
	cp $BrakerGff $FinalDir/final_genes_Braker.gff3
  $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker

  # Combine both fasta files
	cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
	cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
	cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
	cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

  # Combine both gff3 files
	GffBraker=$FinalDir/final_genes_CodingQuary.gff3
	GffQuary=$FinalDir/final_genes_Braker.gff3
	GffAppended=$FinalDir/final_genes_appended.gff3
	cat $GffBraker $GffQuary > $GffAppended

  # Check the final number of genes

	for DirPath in $(ls -d $FinalDir); do
    echo $DirPath;
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
    echo "";
	done
  ```

  ### Remove duplicate and rename genes.

  ```bash
  for GffAppended in $(ls $FinalDir/final_genes_appended.gff3);
  do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    FinalDir=path/to/final/gene/predction/folder #/final
    # Remove duplicated genes
    GffFiltered=$FinalDir/filtered_duplicates.gff
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
    # Rename genes
    GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
    LogFile=$FinalDir/final_genes_appended_renamed.log
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered
    # Create renamed fasta files from each gene feature   
    Assembly=$(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed
    # The proteins fasta file contains * instead of Xs for stop codons, these should be changed
    sed -i 's/\*/X/g' $FinalDir/final/final_genes_appended_renamed.pep.fasta
  done 
```