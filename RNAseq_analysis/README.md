# RNAseq analysis

### Requirements

```bash
conda activate qc_tools
conda install bbmap

conda activate RNAseq
conda install salmon
```

## BBDuk

#### Typical run

```bash
    for RNADir in $(ls -d qc_rna/*/*); do
    FileNum=$(ls $RNADir/F/*_1_trim.fq.gz | wc -l)
        for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
            FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            Ref=/data/scratch/gomeza/prog/bbmap/ribokmers.fa.gz # By courtesy of G.Deakin
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
            echo $StrainPath
            Strain=$(sed 's/.*\///' <<< $StrainPath)
            sbatch -p himem $ProgDir/bbduk.sh $Ref "$StrainPath"/cleaned $FileF $FileR $ProgDir $Strain
        done
    done
```

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

## Transcriptome assembly and differential expression analysis for RNA-Seq with cufflinls


#### Cufflinks. Assembly transcriptomes from RNA-Seq data and quantifies their expression

```bash
    for Alignment in $(ls Fus2_canu_new/Fus2_CzapekDox/*/star_aligmentAligned.sortedByCoord.out.bam); do
        Gff=/projects/fusarium_ex_strawberry/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
        OutDir=$(dirname $Alignment)
        mkdir -p $OutDir/fpkm
        cufflinks -p 8 -o $OutDir/fpkm -G $Gff $Alignment
    done
```

#### Cuffmerge. Merge multiple assembles into one main assembly.

```bash
nano assemblies.txt

    ./Fus2_canu_new/Fus2_CzapekDox/6_S2_L001_R1_001_trim.fq.gz/fpkm/transcripts.gtf
    ./Fus2_canu_new/Fus2_GlucosePeptone/7_S3_L001_R1_001_trim.fq.gz/fpkm/transcripts.gtf
    ./Fus2_canu_new/Fus2_PDA/9_S4_L001_R1_001_trim.fq.gz/fpkm/transcripts.gtf
    ./Fus2_canu_new/Fus2_PDB/4_S1_L001_R1_001_trim.fq.gz/fpkm/transcripts.gtf

    cuffmerge -g /projects/fusarium_ex_strawberry/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 -s /projects/oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa -p 8 assemblies.txt
```

#### Cuffquant. Quantify gene and transcript expression

```bash
    cuffquant -o quant/CzapekDox merged_asm/merged.gtf Fus2_canu_new/Fus2_CzapekDox/*/star_aligmentAligned.sortedByCoord.out.bam
    cuffquant -o quant/GlucosePeptone merged_asm/merged.gtf Fus2_canu_new/Fus2_GlucosePeptone/*/star_aligmentAligned.sortedByCoord.out.bam
    cuffquant -o quant/PDA merged_asm/merged.gtf Fus2_canu_new/Fus2_PDA/*/star_aligmentAligned.sortedByCoord.out.bam
    cuffquant -o quant/PDB merged_asm/merged.gtf Fus2_canu_new/Fus2_PDB/*/star_aligmentAligned.sortedByCoord.out.bam
```

#### Cuffdiff. Compare expression leves of genes and transcripts.

```bash
    cuffdiff merged_asm/merged.gtf Fus2_canu_new/Fus2_CzapekDox/6_S2_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam Fus2_canu_new/Fus2_GlucosePeptone/7_S3_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam Fus2_canu_new/Fus2_PDA/9_S4_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam Fus2_canu_new/Fus2_PDB/4_S1_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam
```

#### Cuffnorm. Normalize expression levesl from a set of RNA-Seq libraries

```bash
    cuffnorm -o norm/v2 --output-format cuffdiff -L Czap,Gluc,PDA,PDB merged_asm/merged.gtf quant/CzapekDox/abundances.cxb quant/GlucosePeptone/abundances.cxb quant/PDA/abundances.cxb quant/PDB/abundances.cxb
```