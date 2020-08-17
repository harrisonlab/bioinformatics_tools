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
    sbatch $ProgDir/count_nucl.sh $F_Read $R_Read 45 #Estimated genome size
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