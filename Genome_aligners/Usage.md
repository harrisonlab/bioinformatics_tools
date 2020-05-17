# Star

## Requirements

```bash
conda install star
```

```bash
for Assembly in $(ls ../assembly/miniasm/N.ditissima/Hg199/racon_10/pilon/filtered_contigs/repeat_masked/Hg199_miniasm_contigs_unmasked.fa)
do
Strain=$(echo $Assembly | rev | cut -f6 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f7 -d '/' | rev)
echo "$Organism - $Strain"
for FileF in $(ls ../qc_rna/RNAseq/N.ditissima/Hg199/mycelium/F/*_trim.fq.gz)
do
FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_trim/_2_trim/g')
echo $FileF
echo $FileR
Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
echo "$Timepoint"
Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
OutDir=alignment_v2/star/$Organism/$Strain
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir
done
done

Alignment outputs were concatenated and Braker1 prediction was run
  Strain=Hg199
  Organism=N.ditissima
  mkdir -p alignment_vAG/star/$Organism/$Strain/mycelium/concatenated
  samtools merge -f alignment_vAG/star/$Organism/$Strain/mycelium/concatenated/concatenated.bam \
  alignment_vAG/star/$Organism/$Strain/mycelium/Hg199_1/star_aligmentAligned.sortedByCoord.out.bam \
  alignment_vAG/star/$Organism/$Strain/mycelium/Hg199_2/star_aligmentAligned.sortedByCoord.out.bam \
  alignment_vAG/star/$Organism/$Strain/mycelium/Hg199_3/star_aligmentAligned.sortedByCoord.out.bam
```