# Amplicon seq analysis

## Potato cyst nematode

### Remove adapters

```bash
# Run fastq-mcf
    for fastq in $(ls ../../../archives/2021_camb_miseq/ANALYSIS/210720_M04543_0009_000000000-DBN4W/Data/Intensities/BaseCalls/*R1_001.fastq.gz); do
        sampleID=$(echo $fastq | rev | cut -f1 -d '/' | rev | sed 's/_S.*//g')
        echo $sampleID
        OutDir=qc_dna/paired/PCN/Plate1/$sampleID
        mkdir -p $OutDir
        Read_F=$fastq
        Read_R=$(echo $Read_F | sed 's/_R1_/_R2_/g')
        echo $Read_F;
        echo $Read_R;
        IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p himem $ProgDir/fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA $OutDir
    done

    for fastq in $(ls ../../../archives/2021_camb_miseq/ANALYSIS/210728_M04543_0011_000000000-DC78R/Data/Intensities/BaseCalls/*R1_001.fastq.gz); do
        sampleID=$(echo $fastq | rev | cut -f1 -d '/' | rev | sed 's/_S.*//g')
        echo $sampleID
        OutDir=qc_dna/paired/PCN/Plate2/$sampleID
        mkdir -p $OutDir
        Read_F=$fastq
        Read_R=$(echo $Read_F | sed 's/_R1_/_R2_/g')
        echo $Read_F;
        echo $Read_R;
        IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p himem $ProgDir/fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA $OutDir
    done
```

### Classification using centrifuge (nt database)

Plate 1

```bash
    for fastq in $(ls ../../../archives/2021_camb_miseq/ANALYSIS/210720_M04543_0009_000000000-DBN4W/Data/Intensities/BaseCalls/*R1_001.fastq.gz); do
    sampleID=$(echo $fastq | rev | cut -f1 -d '/' | rev | sed 's/_S.*//g')
    echo $sampleID
    mkdir -p analysis/Centrifuge/$sampleID
        for Read1 in $(ls ../../../archives/2021_camb_miseq/ANALYSIS/210720_M04543_0009_000000000-DBN4W/Data/Intensities/BaseCalls/"$sampleID"_*R1_001.fastq.gz); do
            Read2=$(echo $Read1 | sed 's/_R1_/_R2_/g')
            echo $Read1
            echo $Read2
            #Sample=$(echo $Read1 | rev | cut -f1 -d '/' | rev | sed 's/_S37_L001_R1_001.fastq.gz//g')
            #echo $Sample
            Database=../../../scratch/public_data/nt-centrifuge_12jan2021/nt
            OutDir=analysis/Centrifuge/$sampleID
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
            sbatch -p short $ProgDir/centrifuge.sh $Database $OutDir $Read1 $Read2
        done
    done
```

Plate 2 

```bash
    for fastq in $(ls ../../../archives/2021_camb_miseq/ANALYSIS/210720_M04543_0009_000000000-DBN4W/Data/Intensities/BaseCalls/*R1_001.fastq.gz); do
    sampleID=$(echo $fastq | rev | cut -f1 -d '/' | rev | sed 's/_S.*//g')
    echo $sampleID
    #mkdir -p analysis/Centrifuge/$sampleID
        for Read1 in $(ls ../../../archives/2021_camb_miseq/ANALYSIS/210720_M04543_0009_000000000-DBN4W/Data/Intensities/BaseCalls/"$sampleID"_*R1_001.fastq.gz); do
            Read2=$(echo $Read1 | sed 's/_R1_/_R2_/g')
            echo $Read1
            echo $Read2
            #Sample=$(echo $Read1 | rev | cut -f1 -d '/' | rev | sed 's/_S37_L001_R1_001.fastq.gz//g')
            #echo $Sample
            Database=../../../scratch/public_data/nt-centrifuge_12jan2021/nt
            OutDir=analysis/Centrifuge/$sampleID
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
            sbatch -p short $ProgDir/centrifuge.sh $Database $OutDir $Read1 $Read2
        done
    done
```

## Slugbot analysis

### Remove adapters

```bash
    for fastq in $(ls ../../../archives/2021_camb_miseq/ANALYSIS/210624_M04543_0008_000000000-DBMTB/Data/Intensities/BaseCalls/*R1_001.fastq.gz); do
        sampleID=$(echo $fastq | rev | cut -f1 -d '/' | rev | sed 's/_S.*//g')
        echo $sampleID
        OutDir=qc_dna/paired/Slugbot/$sampleID
        mkdir -p $OutDir
        Read_F=$fastq
        Read_R=$(echo $Read_F | sed 's/_R1_/_R2_/g')
        echo $Read_F;
        echo $Read_R;
        IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p short $ProgDir/fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA $OutDir
    done
```

### Classification using centrifuge (nt database)


```bash
    for fastq in $(ls qc_dna/paired/Slugbot/*/*/*R1_001_trim.fq.gz); do
        sampleID=$(echo $fastq | rev | cut -f1 -d '/' | rev | sed 's/_S.*//g')
        echo $sampleID
        OutDir=analysis/Centrifuge/Slugbot/$sampleID
        mkdir -p $OutDir
        for Read1 in $(ls qc_dna/paired/Slugbot/*/*/"$sampleID"_*R1_001_trim.fq.gz); do
        Read2=qc_dna/paired/Slugbot/*/R/"$sampleID"_*R2_001_trim.fq.gz
        echo $Read1
        echo $Read2
        #Sample=$(echo $Read1 | rev | cut -f1 -d '/' | rev | sed 's/_S37_L001_R1_001.fastq.gz//g')
        #echo $Sample
        Database=../../../scratch/public_data/nt-centrifuge_12jan2021/nt
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
        sbatch -p long $ProgDir/centrifuge.sh $Database $OutDir $Read1 $Read2
        done
    done

    for kraken in $(ls analysis/Centrifuge/Slugbot/*/centrifuge_krakened.txt); do
        ID=$(echo $kraken | rev | cut -f2 -d '/' | rev )
        echo $ID
        cat $kraken | awk '$4 == "S" { print $0 }' | head -1 > analysis/Centrifuge/Slugbot/"$ID"_stats.txt
    done


    cat analysis/Centrifuge/Slugbot/*_stats.txt > analysis/Centrifuge/Slugbot/Slug_top_hits.txt
```
