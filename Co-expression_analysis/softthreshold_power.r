#!/usr/bin/Rscript

# Load packages

library("WGCNA")
library("optparse")

# Important option recommended in WGCNA documentation

options(stringsAsFactors = FALSE)

# Allow multi-threading

allowWGCNAThreads(nThreads = 4)

# Parse arguments

opt_list <- list(
  make_option("--out_dir", type = "character",
  help = "Directory for plots to be written to"),
  make_option("--max_SFT", type = "integer",
  help = "Maximum SFT value to trial")
  )

opt <- parse_args(OptionParser(option_list = opt_list))
outdir <- opt$out_dir
max_SFT <- opt$max_SFT

# Load previous scripts output

lfile <- paste(outdir, "Cleaned_data.RData", sep = "/")
lnames <- load(file = lfile)

# Step-by-step construction of the network & module ID

# Testing soft-thresholding power values

powers <- c(c(1:max_SFT))

sft <- pickSoftThreshold(datexpr0, powerVector = powers, verbose = 5)

# Draw a plot to allow manual picking of sft value

cex1 <- 0.9
file <- paste(outdir, "sft_testing.pdf", sep = "/")
pdf(file, height = 9, width = 12)
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
xlab = "Soft Threshold (power)",
ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)",
ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1,
col = "red")
dev.off()