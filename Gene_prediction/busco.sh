#!/usr/bin/env bash
#SBATCH -J busco
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=20

### BUSCO analysis to identify single copy genes conserved in Fungi in all genomes in the study. Sample submission script for one genome.
### NOTE: use the genome mode for the transcriptome contigs, otherwise get all contigs re-named as Transcript 1 in the output!

### Do not forget to input the path to the assembly
### And select the DatabaseOpt (options: Eukaryotic, Fungal, Plant, Bacteria)
### NOTE that if you are working with Fungi or Bacteria, you can download lower taxonomic level
### DatabaseOpts to use (e.g. just for Sordariomyceta)

# NOTE to easily prepare a figure with a comparison of complete, fragmented, duplicated
# and missing BUSCO genes in each genome, use the script BUSCO_plot.py in the BUSCO folder.
# Instructions in the user guide BUSCO_v2.0_userguide.pdf
###


Usage="busco.sh [Assembly.fasta] [Eukaryotic/Fungal/Plant/Bacteria or absoulte_path_to_specific_busco_library] [Output_directory]"
echo "$Usage"


Assembly=$1
DatabaseOpt=$2
OutDir=$3

#if [ "$DatabaseOpt" = "Eukaryotic" ]
#then
    #db=/projects/dbBusco/eukaryota_odb9
#elif [ "$DatabaseOpt" = "Fungal" ]
#then
    #db=/projects/dbBusco/fungi_odb9
#elif [ "$DatabaseOpt" = "Plant" ]
#then
    #db=/projects/dbBusco/embryophyta_odb9
#elif [ "$DatabaseOpt" = "Bacteria" ]
#then
    #db=/projects/dbBusco/bacteria_odb9
#elif [ -e $DatabaseOpt ]
#then
    #db=$DatabaseOpt
#else
    #exit
#fi


### Output folder
Filename=$(basename "$Assembly")
Prefix="${Filename%.*}"

### Setting variables
CurDir=$PWD
WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}


### Prep
mkdir -p $WorkDir
cp $Assembly $WorkDir
cd $WorkDir

### Execute
#run_BUSCO.py -o $Prefix -i $Filename -l $db -m geno -c 8 -sp fusarium_graminearum
busco -o $Prefix -i $Filename -l $DatabaseOpt -m geno -c 8 --augustus_species fusarium_graminearum

### Cleanup
rm $Filename
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.
rm -r $WorkDir
