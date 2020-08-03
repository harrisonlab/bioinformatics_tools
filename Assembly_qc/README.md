# Assembly qc tools

Tools used in the quality control and edition of genome assemblies

1. Quast. Quality assessment tool

2. Kat, a K-mer analysis toolkit to assess level of error and duplications in the genome assemblies generated. Mapleson et al., 2016.

3. remove_contaminants.py, split or remove contigs based on tab seperated co-ordinate file and rename contigs.



## Quast

Produce genome assemblies statistics

### Requirements

```bash
# Python 2.7 required
conda install quast
```

### Typical run
```bash
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    for Assembly in $(ls path/to/genome/*.fa); do
        OutDir=$(dirname $Assembly)
        sbatch $ProgDir/quast.sh $Assembly $OutDir
    done
```

## Kat

### Typical run

```bash
    for Assembly in $(ls path/to/genome/assembly/*.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev )
        Organism=$(echo $Assembly | rev | cut -d '/' -f6 | rev)
        echo "$Organism - $Strain"
        IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
        cat $IlluminaDir/F/*_trim.fq.gz > $IlluminaDir/F/F_trim_appended.fq.gz
        cat $IlluminaDir/R/*_trim.fq.gz > $IlluminaDir/R/R_trim_appended.fq.gz
        ReadsF=$(ls $IlluminaDir/F/F_trim_appended.fq.gz)
        ReadsR=$(ls $IlluminaDir/R/R_trim_appended.fq.gz)
        OutDir=$(dirname $Assembly)/kat
        Prefix="${Strain}"
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
        sbatch $ProgDir/kat.sh $Assembly $ReadsF $ReadsR $OutDir $Prefix 200
    done
```

After KAT jobs have finished running, then remove appended trimmed reads

```bash
    rm qc_dna/paired/*/*/*/F_trim_appended.fq.gz
    rm qc_dna/paired/*/*/*/R_trim_appended.fq.gz
```

## Rename contigs

```bash
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    # If split or remove contigs is needed, provide FCSreport file by NCBI.
    touch tmp.txt
    for Assembly in $(ls path/to/assembly/*.fasta); do
        OutDir=$(dirname $Assembly)
        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"GiveGenomeName"_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
    done
    rm tmp.txt
```
