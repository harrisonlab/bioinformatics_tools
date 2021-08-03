# Potato cyst nematode analysis

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




for file in chr*.fa
do
kraken2-build --add-to-library $file --db $DBNAME
done

/scratch/public_data/nt-centrifuge_12jan2021/nt.fa


#### Kraken

```bash
    for Read1 in $(ls raw_data/S3D6-BC11_BC12_No_Selection/barcode1*_filtered.fq.gz); do
        Out=$(echo $Read1 | rev | cut -f1 -d '/' | rev | sed 's/_filtered.fq.gz//g')
        echo $Out
        Database=Foxysporum
        OutDir=analysis/Kraken/filtered_1kb/$Out
        mkdir -p $OutDir
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
        sbatch $ProgDir/kraken_long_reads.sh $Read1 $Database $OutDir
    done

for Reads in $(ls /archives/2021_camb_miseq/ANALYSIS/210720_M04543_0009_000000000-DBN4W/Data/Intensities/BaseCalls); do
Read1=$Reads/M166_S37_L001_R1_001.fastq.gz
Read2=$Reads/M166_S37_L001_R2_001.fastq.gz
centrifuge -x /scratch/public_data/nt-centrifuge_12jan2021/nt -t -1 $Read1 -2 $Read2 --phred33 --report-file centrifuge_report.tsv -S centrifuge_results.txt 
done
```



for fastq in $(ls ../../../archives/2021_camb_miseq/ANALYSIS/210720_M04543_0009_000000000-DBN4W/Data/Intensities/BaseCalls/*R1_001.fastq.gz); do
sample=$(echo $fastq | rev | cut -f1 -d '/' | rev | sed 's/_S.*//g')
echo $sample
mkdir -p analysis/Centrifuge/$sample
done


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

kraken2-build --download-taxonomy --db Nematode
kraken2-build --standard --db Nematode


kraken2-build --standard --db Globodera

kraken2-build --download-library Globodera
centrifuge-download -o library -m -d "Nematode" genbank