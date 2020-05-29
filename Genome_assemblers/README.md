## flye

```bash
  for TrimReads in $(ls qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=WT_minion
    Prefix="$Strain"_flye
    OutDir=assembly/flye/F.venenatum/"$Strain"
    mkdir -p $OutDir
    Size=
    ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/SMARTdenovo
    sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size
  done

  flye --nano-raw qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz --out-dir flye/ --genome-size 37m --threads 8


  RawReads=$1
Prefix=$2
OutDir=$3
Size=$4
```

```bash
# Python 2.7 is needed to install Quast
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc
for Assembly in $(ls flye/assembly.fasta); do
OutDir=$(dirname $Assembly)
sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls flye/assembly.fasta); do
Strain=WT_minion
Organism=F.venenatum
echo "$Organism - $Strain"
ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
done
```
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
racon_round_4.fasta \
qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz \
> racon_round_5.reads_mapped.paf

racon -t 16 qc_dna/minion/F.venenatum/WT/WT_minion_allfiles.fastq.gz racon_round_5.reads_mapped.paf racon_round_4.fasta > racon_round_5.fasta
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

# Build an index mapping from basecalled reads to the signals measured by the sequencer
ScratchDir=/data/scratch/nanopore_tmp_data/$Organism/$Strain
Fast5Dir1=$ScratchDir/path/to/Fast5Files/workspace/pass
Fast5Dir2=$ScratchDir/path/to/Fast5Files/workspace/pass
nanopolish index -d $Fast5Dir1 -d $Fast5Dir2 $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
```

```bash
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_minion_racon_round_10.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ReadDir=../../home/gomeza/rawdata4nanopolish/F.venenatum/WT
OutDir=nanopolish_bwa
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/nanopolish
sbatch $ProgDir/sub_bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
done
```


Split the assembly into 50Kb fragments an submit each to the cluster for nanopolish correction

```bash
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_miniasm_racon10_renamed.fasta ); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly/nanopolish_variants)
RawReads=$(ls ../../home/gomeza/rawdata4nanopolish/F.venenatum/WT_minion/"$Strain"_concatenated_reads_filtered.fastq)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)

#NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
#nanopolish_makerange.py $Assembly > $OutDir/nanopolish/nanopolish_range.txt

Ploidy=1
#echo "nanopolish log:" > nanopolish_log.txt
#for Region in $(cat $OutDir/nanopolish/nanopolish_range.txt | head -n1); do
#Jobs=$(squeue | grep 'nanopo' | grep 'qw' | wc -l)
#while [ $Jobs -gt 1 ]; do
#sleep 1m
#printf "."
#Jobs=$(squeue | grep 'nanopo' | grep 'qw' | wc -l)
#done		
#printf "\n"
#echo $Region
#echo $Region >> nanopolish_log.txt
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/nanopolish
sbatch --array=1-1000%50 $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```



## Rename contigs

Rename the sequences in assembly fasta file to have simple names.

```bash
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.txt
for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_smartdenovo_racon_round_10.fasta); do
OutDir=$(dirname $Assembly)
$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/WT_miniasm_racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
rm tmp.txt
```


Raw reads were moved onto the cluster scratch space for this step and unpacked:
```bash
ScratchDir=/data/scratch/nanopore_tmp_data/Fven
mkdir -p $ScratchDir
cp raw_dna/minion/F.venenatum/WT/*.tar.gz $ScratchDir/.
for Tar in $(ls $ScratchDir/*.tar.gz); do
  tar -zxvf $Tar -C $ScratchDir
done


for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_racon10_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
# if [ -d $ReadDir ]; then
# echo "reads already extracted"
# else
# echo "extracting reads"
mkdir -p $ReadDir
# CurDir=$PWD
# cd $ReadDir
# Event information would have been used from all of the runs, howver MinKnow doesnt
# produce event-level information and therefore just the albacore data was used.
# for Fast5Dir in $(ls -d /home/miseq_data/minion/2017/Fvenenatum/downloaded/pass); do
# nanopolish extract -r $Fast5Dir \
# | gzip -cf
# done > "$Strain"_reads.fa.gz
# cd $CurDir
ReadsFq1=$(ls /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum/raw_dna/minion/F.venenatum/WT/F.venenatum_WT_07-03-17_albacore_v2.02.fastq.gz)
ReadsFq2=$(ls /projects/oldhome/groups/harrisonlab/project_files/fusarium_venenatum/raw_dna/minion/F.venenatum/WT/F.venenatum_WT_18-07-17_albacore_v2.02.fastq.gz)
cat $ReadsFq1 $ReadsFq2 | gunzip -cf > $ReadDir/"$Strain"_concatenated_reads.fastq
echo $ReadsFq1
echo $ReadsFq2

/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_concatenated_reads.fastq --out $ReadDir/"$Strain"_concatenated_reads_filtered.fastq

# cat $ReadsFq1 | gunzip -cf > $ReadDir/"$Strain"_07-03-17_reads.fastq
# cat $ReadsFq2 | gunzip -cf > $ReadDir/"$Strain"_18-07-17_reads.fastq
# /home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_07-03-17_reads.fastq --out $ReadDir/"$Strain"_07-03-17_reads_filtered.fastq
# /home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_18-07-17_reads.fastq --out $ReadDir/"$Strain"_18-07-17_reads_filtered.fastq

ScratchDir=/data/scratch/nanopore_tmp_data/Fven
Fast5Dir1=$ScratchDir/F.venenatum_WT_07-03-17/workspace/pass
Fast5Dir2=$ScratchDir/F.venenatum_WT_18-07-17/workspace/pass
echo $Fast5Dir1
echo $Fast5Dir2
#nanopolish extract -r -q -o $ReadDir/F.venenatum_WT_07-03-17_albacore_v2.02.fastq $Fast5Dir1
nanopolish index -d $Fast5Dir1 -d $Fast5Dir2 $ReadDir/"$Strain"_concatenated_reads_filtered.fastq

# RawReads=$(ls $ReadDir/"$Strain"_reads_no_duplicates.fa.gz)
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
# submit alignments for nanoppolish
sbatch $ProgDir/bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
done





for Assembly in $(ls ../fusarium_ex_pea/assembly/SMARTdenovo/*/*/racon/racon_min_500bp_renamed.fasta | grep 'FOP1'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq1=$(ls ../fusarium_ex_pea/raw_dna/minion/*/FOP1-EMR/*.fastq.gz)
cp $ReadsFq1 $ReadDir/.
Reads=$(ls $ReadDir/*.fastq.gz)

ScratchDir=/data/scratch/nanopore_tmp_data/Fven
Fast5Dir1=$ScratchDir/F.oxysporum_fsp_pisi_FOP1-EMR_2017-07-11/workspace/pass
nanopolish index -d $Fast5Dir1 $ReadsFq1
# OutDir=$(dirname $Assembly)
OutDir=$ReadDir
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $Reads $OutDir/nanopolish
done







Split the assembly into 50Kb fragments an submit each to the cluster for nanopolish correction

for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_minion_racon10_renamed.fasta); do
Assembly=assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_minion_racon10_renamed.fasta
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)
RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/"$Strain"_concatenated_reads_filtered.fastq)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)
Assembly=assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_minion_racon10_renamed.fasta
#NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
nanopolish_makerange.py $Assembly > $OutDir/nanopolish/nanopolish_range.txt
Ploidy=1

un this one!!!!


i=1
for Region in $(cat $OutDir/nanopolish/nanopolish_range.txt)
do
cat $Region | sed -n ''$i'p'
printf "\n"
echo $Region
echo $Region >> nanopolish_log.txt
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
$ProgDir/nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/region/$Region
i=$(($i+1))
done

nanopolish_makerange.py assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_minion_racon10_renamed.fasta | parallel --results nanopolish.results -P 8 \
nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r raw_dna/nanopolish/F.venenatum/WT_minion/WT_minion_concatenated_reads_filtered.fastq \
--max-haplotypes 100000 \
--fix-homopolymers \
--min-candidate-frequency 0.2 \
-b assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/nanopolish/reads.sorted.bam \
-g assembly/SMARTdenovo/F.venenatum/WT_minion/racon_10/WT_minion_racon10_renamed.fasta -t 4

nanopolish variants \
  -t 8 \
  --ploidy 1 \
  -w contig_1 \
  --consensus \
  --max-haplotypes 100000 \
  --fix-homopolymers \
  --min-candidate-frequency 0.2 \
  --reads raw_dna/nanopolish/F.venenatum/WT_minion/WT_minion_concatenated_reads_filtered.fastq \
  --bam assembly/miniasm/F.venenatum/WT_minion/racon_10/nanopolish/reads.sorted.bam \
  --genome assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_miniasm_racon10_renamed.fasta \
  > consensus_variants.txt




i=$(($i+1))
done

echo "nanopolish log:" > nanopolish_log.txt
for Region in $(cat $OutDir/nanopolish/nanopolish_range.txt | sed -n ''$i'p')
do
printf "\n"
echo $Region
i=$(($i+1))
done
done
done
echo $Region >> nanopolish_log.txt
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
$ProgDir/nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
i=$(($i+1))
done

i=1
for FILE in $(cat assembly/miniasm/F.venenatum/WT_minion/racon_10/nanopolish/nanopolish_range.txt)
do
sed -n ''$i'p'
printf "\n"
i=$(($i+1))


#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
#$ProgDir/nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
i=$(($i+1))
done
>> nanopolish_log.txt
srun --partition medium -mem 10G --cpus-per-task 10 --pty bash

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

for Assembly in $(ls assembly/SMARTdenovo/F.venenatum/WT_albacore_v2/nanopolish/WT_albacore_v2_nanoplish_min_500bp_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
	# Quast
  OutDir=$(dirname $Assembly)
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
	# Busco
	BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
	OutDir=gene_pred/busco/$Organism/$Strain/assembly
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
	qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done