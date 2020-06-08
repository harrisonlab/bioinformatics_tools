# Genome assemblers

Genome assemblers

1. 

2. 

3. 

4. 



## Flye

De novo assembler for single molecule sequencing reads using repeat graphs

### Requirements

```bash
conda install flye
```

### Typical run

```bash
  for TrimReads in $(ls path/to/single/molecule/sequencing/reads/*_allfiles.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
    Prefix="$Strain"_flye
    OutDir=assembly/flye/$Organism/$Strain
    mkdir -p $OutDir
    Size=37m # Expected genome size
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools
    sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size
  done
```
  
## Racon


```bash
for Assembly in $(ls flye/assembly.fasta); do
ReadsFq=$(ls qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

minimap2 \
-x map-ont \
-t16 \
racon_round_9.fasta \
qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz \
> racon_round_10.reads_mapped.paf

racon -t 16 qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz racon_round_10.reads_mapped.paf racon_round_9.fasta > racon_round_10.fasta
cp racon_round_$i.fasta current-assembly.fa
cp racon_round_$i.fasta $CurDir/$OutDir/"$Prefix"_racon_round_$i.fasta


Assembly correction using nanopolish
Fast5 files are very large and need to be stored as gzipped tarballs. These needed temporarily unpacking but must be deleted after nanpolish has finished running.


faidx -d '|' final_genes_appended_renamed.cdna.fasta $(tr '\n' ' ' < Tri5_genelist.txt) > selected_genes.fasta



## Assembly correction with nanopolish

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
    OutDir=nanopolish_bwa # Output directory
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






```bash



for Assembly in $(ls assembly/SMARTdenovo/*/*/racon/racon_min_500bp_renamed.fasta | grep 'WT' | grep 'albacore'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
InDir=$(dirname $Assembly)
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_merge.py $InDir/*:*-*/*.fa > $OutDir/"$Strain"_nanoplish.fa

echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
Quast and busco were run to assess the effects of nanopolish on assembly quality:



```bash
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.txt
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_minion_racon_round_10.fasta); do
  OutDir=$(dirname $Assembly)
  $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/WT_miniasm_racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
rm tmp.txt
```

Medaka


conda create -n medaka -c conda-forge -c bioconda medaka


```bash
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/WT_minion_miniasm.fa); do
ReadsFq=$(ls qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz)
OutDir=$(dirname $Assembly)/medaka3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
done
```
medaka_consensus -i qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz -d assembly/miniasm/F.venenatum/WT_minion/WT_minion_miniasm.fa -o medaka2 -t 8

```bash
for Assembly in $(ls assembly/flye/F.venenatum/WT_minion/WT_minion_miniasm.fa); do
ReadsFq=$(ls qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz)
OutDir=$(dirname $Assembly)/medaka3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
done
```