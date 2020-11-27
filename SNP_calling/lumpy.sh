#!/usr/bin/env bash
#SBATCH -J lumpy
#SBATCH --partition=long 
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=10


Prefix=$1
Assembly=$2
F_reads=$3
R_reads=$4
OutDir=$5

CWD=$PWD

WorkDir="$TMPDIR"
mkdir -p $WorkDir

cp -r $Assembly $F_reads $R_reads $WorkDir
Fr=$(basename "$F_reads")
Rr=$(basename "$R_reads")
As=$(basename "$Assembly")

cd $WorkDir
Output=${Fr%.*}.bam

# Align the data
bwa index $As

bwa mem -t 12 $As $Fr $Rr | samtools view -S -b - > "$Prefix".bam


bamaddrg=/projects/oldhome/sobczm/bin/freebayes/bamaddrg/bamaddrg
Output_rg_un=${Output%.bam}_rg_unsoerted.bam
Output_rg=${Output%.bam}_rg.bam
### Add group and sample name (Prefix)
$bamaddrg -b "$Prefix".bam -s $Prefix -r $Prefix > "$Prefix"_unsorted.bam
### Sort the full BAM file.
samtools sort "$Prefix"_unsorted.bam -o "$Prefix"_sorted.bam

# Extract the discordant paired-end alignments.
samtools view -b -F 1294 "$Prefix"_sorted.bam > "$Prefix"_discordants_unsorted.bam

# Extract the split-read alignments
samtools view -h "$Prefix"_sorted.bam | extractSplitReads_BwaMem -i stdin | samtools view -Sb - > "$Prefix"_splitters_unsorted.bam

# Sort both alignments
samtools sort "$Prefix"_discordants_unsorted.bam -o "$Prefix"_discordants_sorted.bam
samtools sort "$Prefix"_splitters_unsorted.bam -o "$Prefix"_splitters_sorted.bam

# printf "@HD\tVN:1.5\tSO:coordinate\n" > $TMPDIR/sorted_line.txt
# for BAM in $Output_rg ${Output_rg%.bam}_discordants.bam ${Output_rg%.bam}_splitters.bam
# do
#      samtools view -H $BAM > header.sam
#      cat header.sam $TMPDIR/sorted_line.txt>new_header.sam
#      samtools reheader new_header.sam $BAM >${BAM%.bam}_test.bam
#      mv ${BAM%.bam}_test.bam $BAM
# done


#index
samtools index "$Prefix"_sorted.bam
samtools index "$Prefix"_discordants_sorted.bam
samtools index "$Prefix"_splitters_sorted.bam

rm ${Output_rg%.bam}_splitters_unsorted.bam ${Output_rg%.bam}_discordants_unsorted.bam $Output $Output_rg_un
rm $Fr $Rr
mkdir -p $CWD/$OutDir
cp -r *_sorted* $CWD/$OutDir/.
rm -rf $WorkDir
