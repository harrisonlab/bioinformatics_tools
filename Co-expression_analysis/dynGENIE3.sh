#!/usr/bin/env bash
#SBATCH -J dynGenie3
#SBATCH --partition=long
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4

scripts=/mnt/shared/scratch/agomez/apps/git_repos/bioinformatics_tools/Co-expression_analysis

/mnt/shared/scratch/agomez/apps/conda/envs/Rstudio/bin/Rscript --vanilla $scripts/dynGENIE3_script.r --exp1 $1 --exp2 $2 --exp3 $3 --exp4 $4 --exp5 $5 --tfs $6 --out_dir $7