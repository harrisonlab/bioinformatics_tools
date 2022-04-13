#!/mnt/shared/scratch/agomez/apps/conda/envs/perly_env/bin/Rscript

# Load libraries

library(optparse)
library(reshape2)
library(doRNG)
library(doParallel)
library(dplyr)
source ("dynGENIE3_R_C_wrapper/dynGENIE3.R")

# Parse in arguments to script

opt_list <- list(
    make_option("--exp1", type = "character",
    help = "vst experiment 1"),
    make_option("--exp2", type = "character",
    help = "vst experiment 2"),
    make_option("--exp3", type = "character",
    help = "vst experiment 3"),
    make_option("--exp4", type = "character",
    help = "vst experiment 4"),
    make_option("--exp5", type = "character",
    help = "vst experiment 5"),
    make_option("--tfs", type = "character",
    help = "regulators list"),
    make_option("--out_dir", type = "character",
    help = "Directory for output to be written to")
  )

opt <- parse_args(OptionParser(option_list = opt_list))
T1 <- opt$exp1
T2 <- opt$exp2
T3 <- opt$exp3
T4 <- opt$exp4
T5 <- opt$exp5
TF <- opt$tfs
outdir <- opt$out_dir

# Edit headers before load matrices into R with read.expr.matrix function (from dynGenie3)
TS1 <- read.expr.matrix(T1,form="rows.are.genes")
TS2 <- read.expr.matrix(T2,form="rows.are.genes")
TS3 <- read.expr.matrix(T3,form="rows.are.genes")
TS4 <- read.expr.matrix(T4,form="rows.are.genes")
TS5 <- read.expr.matrix(T5,form="rows.are.genes")

# Time
time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,], TS5[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),], TS5[2:nrow(TS5),])

# Add regulators
TFexp<-read.table(TF,header=TRUE,sep="\t")
TF<-TFexp[,1]

# Default but needed
decay<-0.02
# Add to reproduce results
set.seed(123)
# Use the Extra-Trees as tree-based method
tree.method <- "ET"
# Number of randomly chosen candidate regulators at each node of a tree
K <- "all"
# Number of trees per ensemble
ntrees <- 500
# Run the method with these settings
resall <- dynGENIE3(TS.data,time.points, regulators=TF, tree.method=tree.method, K=K, ntrees=ntrees,alpha=decay)