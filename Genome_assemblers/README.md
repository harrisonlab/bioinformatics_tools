# Genome assemblers

Genome assemblers and correction tools.

1. Flye: de novo assembler for single molecule sequencing read using repeat graphs.

2. SMARTdenovo: OLC-based de novo assembler for uncorrected long reads using a single perl script.

3. Miniasm: OLC-based de novo assembler for uncorrected long reads. Minimap2 is needed to perform the all-vs-all read self-mapping.

4. Racon: Correct raw contigs from rapid assembly methods, such as Flye, SMARTdenovo and Miniasm.

5. Medaka: Quick and accurate assembly correction using graph-based method.

6. Nanopolish: Correction of assembled nanopore reads, detection of base modification and sequence variants with respect to a reference genome. 

7. Pilon: Polish genome assemblies using read aligments to identify inconsistencies between reads and input genome.

7. Canu


Quast, BUSCO and Kat can be used to evaluate the quality of the assemblies.



### Requirements

```bash
conda activate olc_assemblers
# Flye
conda install flye
# SMARTdenovo
conda install smartdenovo
# Miniasm
conda install miniasm
conda install bbmap # Rename reads with unique ID
conda install minimap2 # all-vs-all mappings
# Racon
conda install racon
# Medaka. This will create an env for medaka only.
conda create -n medaka -c conda-forge -c bioconda medaka
# Nanopolish. Medaka is recommended.
conda install nanopolish 
# Pilon
conda install pilon
```

## Typical run

### Flye

```bash
  for TrimReads in $(ls path/to/single/molecule/sequencing/reads/*_allfiles.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev) # Edit to set your ouput directory
    Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev) # Edit to set your ouput directory
    Prefix="$Strain"_flye
    OutDir=assembly/flye/$Organism/$Strain
    mkdir -p $OutDir
    Size=37m # Expected genome size
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools
    sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size
  done
```
  
### SMARTdenovo

```bash
  for TrimReads in $(ls path/to/single/molecule/sequencing/reads/*_allfiles.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev) # Edit to set your ouput directory
    Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev) # Edit to set your ouput directory
    Prefix="$Strain"_smartdenovo
    OutDir=assembly/SMARTdenovo/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools
    sbatch $ProgDir/SMARTdenovo.sh $TrimReads $Prefix $OutDir
  done
```

### Miniasm

```bash
  for TrimReads in $(ls path/to/single/molecule/sequencing/reads/*_allfiles.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev) # Edit to set your ouput directory
    Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev) # Edit to set your ouput directory
    Prefix="$Strain"_miniasm
    OutDir=assembly/miniasm/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools
    sbatch $ProgDir/miniasm.sh $TrimReads $Prefix $OutDir
  done
```

## Racon

Consensus module for raw de novo DNA assembly

```bash
  for Assembly in $(ls path/to/raw/assembled/contigs/*.fasta); do
    ReadsFq=$(ls path/to/single/molecule/sequencing/reads/*_allfiles.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
# You might rename your contigs at this point using remove_contaminants.py
```

## Medaka


```bash
  for Assembly in $(ls path/to/corrected/consensus/reads/racon_10/*racon10_renamed.fasta); do
    ReadsFq=$(ls path/to/single/molecule/sequencing/reads/*_allfiles.fastq.gz)
    OutDir=$(dirname $Assembly)/medaka
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
  done
```

## Nanopolish

Since Medaka is recommended over Nanopolish for assembly correction, the next commands might not be necessary. 


```bash
  ReadDir=rawdata4nanopolish/$Organism/$Strain
  mkdir -p $ReadDir
  ReadsFq1=$(ls /path/to/raw/basecalled/minion/reads/e.g./F.venenatum_WT_07-03-17_albacore_v2.02.fastq.gz)
  ReadsFq2=$(ls /path/to/raw/basecalled/minion/reads/e.g./F.venenatum_WT_18-07-17_albacore_v2.02.fastq.gz)
  cat $ReadsFq1 $ReadsFq2 | gunzip -cf > $ReadDir/"$Strain"_concatenated_reads.fastq

# Remove duplicate reads

  /home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assembler/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_concatenated_reads.fastq --out $ReadDir/"$Strain"_concatenated_reads_filtered.fastq

# Raw reads were moved onto the cluster scratch space for this step and unpacked

  ScratchDir=/data/scratch/nanopore_tmp_data/$Organism/$Strain
  mkdir -p $ScratchDir
  cp path/to/raw/basecalled/reads/*.tar.gz $ScratchDir/.
  for Tar in $(ls $ScratchDir/*.tar.gz); do
    tar -zxvf $Tar -C $ScratchDir
  done

# Build an index mapping from basecalled reads to the signals measured by the sequencer

  ScratchDir=/data/scratch/nanopore_tmp_data/$Organism/$Strain
  Fast5Dir1=$ScratchDir/path/to/Fast5Files/workspace/pass
  Fast5Dir2=$ScratchDir/path/to/Fast5Files/workspace/pass
  nanopolish index -d $Fast5Dir1 -d $Fast5Dir2 $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
```

### Alignment of minion reads to a minion assembly

```bash
  for Assembly in $(ls path/to/consensus/assembly/*_racon_round_10_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ReadDir=rawdata4nanopolish/$Organism/$Strain # Path to fastq reads
    OutDir=$(dirname $Assembly)/nanopolish_bwa # Output directory
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/nanopolish
    sbatch $ProgDir/bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
  done
```

### Detect SNPs and indels with respect to a reference genome and calculate an improved consensus sequence

```bash
  for Assembly in $(ls path/to/consensus/assembly/*_racon_round_10_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Assembly)/nanopolish_bwa
    RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/"$Strain"_concatenated_reads_filtered.fastq)
    AlignedReads=$(ls $OutDir/reads.sorted.bam)
    Ploidy=1
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $OutDir/variants
  done
```


### Merge nanopolish results 

```bash
  for Assembly in $(ls path/to/consensus/assembly/*_racon_round_10_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    OutDir=$(dirname $Assembly)/nanopolish_final
    mkdir -p $OutDir
    InDir=$(dirname $Assembly)
    NanoPolishDir=/home/gomeza/miniconda3/envs/olc_assemblers/bin # Path to your conda installation path
    python $NanoPolishDir/nanopolish_merge.py $InDir/*:*-*/*.fa > $OutDir/"$Strain"_nanoplish.fa
    echo "" > tmp.txt
    # Rename contigs
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    $ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
  done
```

## Pilon 


```bash
  for Assembly in $(ls path/to/genome/assembly/*_renamed.fasta); do
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  IlluminaDir=$(ls -d path/to/qc_dna/paired/$Organism/$Strain)
  echo $Strain
  echo $Organism
  TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
  TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
  TrimF2_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
  TrimR2_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
  echo $TrimF1_Read
  echo $TrimR1_Read
  echo $TrimF2_Read
  echo $TrimR2_Read
  OutDir=$(dirname $Assembly)
  Iterations=10
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon # Scripts for 1,2 or 3 libraries
  sbatch $ProgDir/sub_pilon_2_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir $Iterations
  done
  # You might rename your contigs at this point using remove_contaminants.py
```

