#!/usr/bin/Rscript

# Load libraries

library("WGCNA")
library("optparse")

# Important option recommended by WGCNA

options(stringsAsFactors = FALSE)

# Parse in arguments to script

opt_list <- list(
  make_option("--gene_table", type = "character",
  help = "Input file of RNA-Seq data"),
  make_option("--out_dir", type = "character",
  help = "Directory for plots to be written to"),
  make_option("--column_start", type = "integer",
  help = "Column number for the start of normalized values"),
  make_option("--column_end", type = "integer",
  help = "Column number for the end of normalized values")
  )

opt <- parse_args(OptionParser(option_list = opt_list))
inp <- opt$gene_table
outdir <- opt$out_dir
column_start <- opt$column_start
column_end <- opt$column_end

# Load input file

exp_data <- read.csv(inp, sep = "\t")

# Parse data as WGCNA tutorial recommends

datexpr0 <- as.data.frame(t(exp_data[, c(column_start:column_end)]))
names(datexpr0) <- exp_data$ID # trancript_ID column
rownames(datexpr0) <- names(exp_data)[c(column_start:column_end)]

# Check for excessive missing values and ID outliers

gsg <- goodSamplesGenes(datexpr0, verbose = 3)
gsg$allOK

# Remove any genes and samples that do not pass the cut

if (!gsg$allOK){
    # Print items removed to file
    if (sum(!gsg$goodGenes) > 0){
        genes_to_remove <- (paste("Removing genes:", paste(names(datexpr0)[
        !gsg$goodGenes], collapse = "\n")))
        gfile <- paste(outdir, "removed_genes.txt", sep = "/")
        write(genes_to_remove, file = gfile)
    }
    if (sum(!gsg$goodSamples) > 0){
        samples_to_remove <- (paste("Removing samples:",
        paste(rownames(datexpr0)[!gsg$goodSamples], collapse = "\n")))
        sfile <- paste(outdir, "removed_samples.txt", sep = "/")
        write(samples_to_remove, file = sfile)
    }
    # Remove items that fail QC
    datexpr0 <- datexpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster samples to check for outliers

sampletree <- hclust(dist(datexpr0), method = "average")
file <- paste(outdir, "sample_clustering.pdf", sep = "/")
pdf(file, width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampletree, main = "Sample clustering to detect outliers", sub = "",
xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Remove outlier samples, the height may need changing so be sure to check

abline(h = 30000, col = "red")
clust <- cutreeStatic(sampletree, cutHeight = 30000, minSize = 10)
table(clust)
keepsamples <- (clust == 0)
datexpr <- datexpr0[keepsamples, ]
ngenes <- ncol(datexpr)
nsamples <- nrow(datexpr)
dev.off()

Rfile <- paste(outdir, "Cleaned_data.RData", sep = "/")
save(datexpr, file = Rfile)
