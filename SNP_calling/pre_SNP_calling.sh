#!/usr/bin/env bash
#SBATCH -J pre_snp_calling
#SBATCH --partition=long
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16


##############################################
# Prep mappings from Bowtie2 for SNP calling
### Remove multimapping reads, discordant reads. PCR and optical duplicates, and
### add read group and sample name to each mapped read (preferably, the shortest ID possible)
#INPUT:
# 1st argument: input SAM file with your mappings
# 2nd argument: sample name (prefix) to be used to identify it in the future
#OUTPUT:
# Indexed BAM file with suffix "nodup_rg" to be fed into SNP calling with GATK.
#############################################
input_sam=$1
prefix=$2
OutDir=$3
filename=$(basename "$input_sam")
name="${filename%.*}"

CurPath=$PWD
WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}

### Prep

WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cp $input_sam $WorkDir
cd $WorkDir

### Get rid of multimapping reads by filtering out on the XS:i: tag
grep -v "XS:i" $filename >temp && mv temp $filename
samtools view -bS -o $name.bam $filename
samtools sort $name.bam $name\_sorted
samtools index $name\_sorted.bam

### Keep only reads with "paired reads" and "properly paired reads" flags.
samtools view -b -h -f 3 -o $name\_proper.bam $name\_sorted.bam
### Sort for downstream analyses
samtools sort $name\_proper.bam $name\_proper\_sorted
samtools index $name\_proper\_sorted.bam

### Remove PCR and optical duplicates
#picard=/home/sobczm/bin/picard-tools-2.5.0/picard.jar
#java -jar $picard MarkDuplicates \
picard MarkDuplicates \
INPUT=$name\_proper\_sorted.bam \
OUTPUT=$name\_proper\_sorted\_nodup.bam \
METRICS_FILE=$name\_proper\_sorted\_nodup.txt \
REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE MAX_RECORDS_IN_RAM=500000000 \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT
### Add group and sample name (prefix)
#java -jar $picard AddOrReplaceReadGroups \
picard AddOrReplaceReadGroups \
INPUT=$name\_proper\_sorted\_nodup.bam \
OUTPUT=$name\_proper\_sorted\_nodup_rg.bam \
SORT_ORDER=coordinate CREATE_INDEX=true RGID=$prefix  RGSM=$prefix \
RGPL=Illumina RGLB=library RGPU=barcode VALIDATION_STRINGENCY=LENIENT
samtools index $name\_proper\_sorted\_nodup_rg.bam

### Cleanup
mv $name\_proper\_sorted\_nodup.txt $CurPath/$OutDir/.
mv $name\_proper\_sorted\_nodup_rg.bam $CurPath/$OutDir/.
mv $name\_proper\_sorted\_nodup_rg.bam.bai $CurPath/$OutDir/.
rm -rf $WorkDir
