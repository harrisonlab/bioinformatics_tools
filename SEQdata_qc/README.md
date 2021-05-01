# QC tools for multi-platform sequence data

1. FastQC: Quality control checks on raw sequence data

2. fast-mcf: Remove adapters and low quality reads

3. Porechop: Adapter removal in Oxford Nanopore reads


### Requirements

```bash
conda activate qc_tools
conda install fastqc
conda install ea-utils
# Porechop is installed in /scratch/software/. Add this line to your profile.
PATH=${PATH}:/scratch/software/Porechop-0.2.3
. ~/.profile # Refresh profile
```

### Typical run

```bash
# Run fastqc
for Strain in Strain1 Strain2; do
    RawData=$(ls raw_dna/paired/$Organism/$Strain/*/*.fastq.gz)
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/fastqc.sh $RawData
done
```

```bash
# Run fastq-mcf
for Strain in Strain1 Strain2; do
    Read_F=raw_dna/paired/*/*/F/*.fastq.gz
    Read_R=raw_dna/paired/*/*/R/*.fastq.gz
    echo $Read_F;
    echo $Read_R;
    IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

```bash
# Estimate coverage (optional)
for DataDir in $(ls -d raw_dna/paired/$Organism/$Strain); do
    F_Read=$(ls $DataDir/F/*.gz)
    R_Read=$(ls $DataDir/R/*.gz)
    echo $F_Read
    echo $R_Read
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/count_nucl.sh $F_Read $R_Read 45 $DataDir #Estimated genome size and DataDir as output directory
done

# Estimate coverage long read data
for RawData in $(ls -d raw_dna/minion/$Organism/$Strain/*fq.gz); do
    echo $RawData
    GenomeSize=45 #Estimated genome size
    OutDir=$(dirname $RawData)
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/count_nucl_single.sh $RawData $GenomeSize $OutDir
done
```


```bash
# Adapter removal ONT reads
    RawReads=path/to/ONT/raw/reads/*.fastq.gz
    OutDir=qc_dna/minion/$Organism/$Strain #e.g.
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/porechop.sh $RawReads $OutDir 
```

## Other commands

```bash
    # Count reads in fastq files
    echo $(zcat file.fq.gz|wc -l)/4|bc
    # Length per read
    python fastq_lenght.py file.fastq 
    # Read length occurrence
    awk '{if(NR%4==2) print length($1)}' file.fastq | sort -n | uniq -c > read_length.txt
    # or 
    zcat file.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > read_length.txt
```
```r
    # Plot read occurrence
    reads<-read.csv(file="read_length.txt", sep="", header=FALSE)
    plot (reads$V2,reads$V1,type="l",xlab="read length",ylab="occurences",col="blue")
```


## Coverage analysis

```bash
  for Bam in $(ls alignment/bwa/vs_*/$Organism/$Strain/*sorted.bam); do
    Strain=$(echo $Bam | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Bam)
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_*_depth.tsv
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_*_depth.tsv > $OutDir/${Organism}_${Strain}_vs_*_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_*_depth_10kb.tsv
  done
  OutDir=analysis/genome_alignment/bwa/vs_*_Ref/grouped
  mkdir -p $OutDir
  cat alignment/bwa/vs_*/$Organism/$Strain/*_*_vs_*_depth_10kb.tsv > analysis/genome_alignment/bwa/vs_*_Ref/grouped/vs_*_grouped_depth.tsv
```