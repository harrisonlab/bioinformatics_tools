# RNAseq analysis

## Salmon 

#### Typical run

```bash
    for Transcriptome in $(ls path/to/predicted/transcriptome/final_genes_appended_renamed.cdna.fasta); do
        Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
        Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
        echo "$Organism - $Strain"
        for RNADir in $(ls -d path/to/RNAseq/reads); do
        FileF=$(ls $RNADir/*.1.fq) # grep -e 'Sample name'
        FileR=$(ls $RNADir/*.2.fq) # grep -e 'Sample name'
        echo $FileF
        echo $FileR
        Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/.1.fq//g')
        echo $Sample_Name
        OutDir=RNAseq_analysis/salmon/$Organism/$Strain/$Sample_Name
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
        sbatch $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
        done
    done
```