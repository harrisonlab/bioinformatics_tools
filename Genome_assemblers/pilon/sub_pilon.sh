#!/usr/bin/env bash
#SBATCH -J pilon
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=30

# Align raw reads to a pacbio assembly and then use this alignmeant to correct
# indels and substitutions in the assembly.

Mem="470G"
#Mem="70G"
# Mem="372G"
Threads=8

# ---------------
# Step 1
# Collect inputs
# ---------------

Assembly=$(basename $1)
Read_F=$(basename $2)
Read_R=$(basename $3)
OutDir=$4
Iterations=$5
if [ $6 ]; then
  Ploidy=$6
else
  Ploidy="haploid"
fi

CurDir=$PWD
echo  "Running Pilon with the following inputs:"
echo "Pacbio assembly - $Assembly"
echo "Forward trimmed reads - $Read_F"
echo "Reverse trimmed reads - $Read_R"
echo "OutDir - $OutDir"
echo "Running Pilon the following number of times - $Iterations"
echo "Ploidy set to: $Ploidy"

mkdir -p $CurDir/$OutDir

# ---------------
# Step 2
# Copy data
# ---------------

WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir
cp $CurDir/$1 assembly.fa
cp $CurDir/$2 $Read_F
cp $CurDir/$3 $Read_R

mkdir best_assembly
cp assembly.fa best_assembly/.

for i in $(seq 1 $Iterations); do
  echo "Running Iteration: $i"
  mkdir $WorkDir/"correction_$i"
  cd $WorkDir/correction_$i
  cp $WorkDir/best_assembly/assembly.fa .

  # ---------------
  # Step 3.a
  # Align seq reads
  # ---------------
  # Prepare the assembly for alignment
  # Align reads against the assembly
  # Convert the SAM file to BAM in preparation for sorting.
  # Sort the BAM file, in preparation for SNP calling:
  # Index the bam file


  bowtie2-build assembly.fa assembly.fa.indexed
  bowtie2 -p 12 -x assembly.fa.indexed -1 $WorkDir/$Read_F -2 $WorkDir/$Read_R  -S assembly.fa_aligned.sam
  samtools view --threads 20 -bS assembly.fa_aligned.sam -o assembly.fa_aligned.bam
  samtools sort --threads 20 assembly.fa_aligned.bam -o assembly.fa_aligned_sorted.bam
  samtools index assembly.fa_aligned_sorted.bam

  # ---------------
  # Step 3.b
  # Run Pilon
  # ---------------
  # Run pilon to polish
  if [ $Ploidy == "haploid" ]; then
    JavaDir=/scratch/software/pilon-1.23
    java -Xmx$Mem -jar $JavaDir/pilon-1.23.jar --threads 16 --genome assembly.fa --changes --frags assembly.fa_aligned_sorted.bam --outdir .
  elif [ $Ploidy == "diploid" ]; then
     JavaDir=/scratch/software/pilon-1.23
    java -Xmx$Mem -jar $JavaDir/pilon-1.23.jar --threads 16 --genome assembly.fa --changes --diploid --frags assembly.fa_aligned_sorted.bam --outdir .
  else
  echo "ploidy not recognised"
  fi
  cp pilon.fasta $WorkDir/best_assembly/assembly.fa
  # cp pilon.changes $WorkDir/best_assembly/pilon_$i.changes
  cp pilon.fasta $CurDir/$OutDir/pilon_$i.fasta
  cp pilon.changes $CurDir/$OutDir/pilon_$i.changes
  cd $WorkDir
done

mv $WorkDir/best_assembly/assembly.fa $WorkDir/best_assembly/pilon.fasta

# mkdir -p $CurDir/$OutDir
# cp $WorkDir/best_assembly/* $CurDir/$OutDir/.
rm -r $WorkDir
