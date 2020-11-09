for Assembly in $(ls race_1_smartdenovo_racon_round_10_renamed.fasta); do
BAM=akin/bwa_mapping.SDEN.sorted.bam
OutDir=akin
Iterations=2
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
sbatch $ProgDir/sub_pilon.sh $Assembly $BAM $OutDir $Iterations
done

srun --partition himem --mem-per-cpu 20G --cpus-per-task 40 --pty bash


java -Xmx378G -jar $JavaDir/pilon-1.17.jar --threads 16 --genome assembly.fa --changes --bam assembly.fa_aligned_sorted.bam --outdir .

```bash
for Assembly in $(ls race_1_smartdenovo_racon_round_10_renamed.fasta); do
Organism=F.oxysporum_fsp_lactucae
Strain=AJ520
IlluminaDir=$(ls -d $Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
echo $TrimF1_Read
echo $TrimR1_Read
OutDir=$(dirname $Assembly)
Iterations=10
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon 
sbatch $ProgDir/sub_pilon_v2.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
  # You might rename your contigs at this point using remove_contaminants.py
```

### BUSCO
```bash
  for Assembly in $(ls pilon_10.fasta); do
    Organism=F.oxysporum_fsp_lactucae
    Strain=AJ520
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10/pilon_10
    sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
  done
```

