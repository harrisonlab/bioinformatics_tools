#!/usr/bin/Rscript

# Load libraries

library("WGCNA")
library("optparse")

# Important option recommended in WGCNA documentation

options(stringsAsFactors = FALSE)

# Allow multi-threading

allowWGCNAThreads(nThreads = 4)

# Parse arguments

opt_list <- list(
    make_option("--out_dir", type = "character",
    help = "Directory for output to be written to"),
    make_option("--sft", type = "integer",
    help = "Value of sft identified from choose_softthreshold.R"),
    make_option("--min_module_size", type = "integer",
    help = "Minimum module size for cutting clustered tree"),
    make_option("--merging_threshold", type = "double",
    help = "Threshold to merge modules with correlation of 1 minus
    specified value.")
    )

opt <- parse_args(OptionParser(option_list = opt_list))
outdir <- opt$out_dir
softpower <- opt$sft
min_mod_size <- opt$min_module_size
medissthres <- opt$merging_threshold

lfile <- paste(outdir, "Cleaned_data.RData", sep = "/")
lnames <- load(file = lfile)

# Calculate adjacency

adjacency <- adjacency(datexpr0, power = softpower)

# Topological Overlap Matrix (TOM)

tom <- TOMsimilarity(adjacency)
disstom <- 1 - tom

# Clustering using TOM

genetree <- hclust(as.dist(disstom), method = "average")

file <- paste(outdir, "clustering_tree.pdf", sep = "/")
pdf(file, height = 9, width = 12)
sizeGrWindow(12, 9)
plot(genetree, xlab = "", sub = "", main = "Gene clustering on TOM-based
dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

# Cut clustering tree into several modules

dynamicmods <- cutreeDynamic(dendro = genetree, distM = disstom, deepSplit = 2,
pamRespectsDendro = FALSE, minClusterSize = min_mod_size)
table(dynamicmods)

# Plot modules on clustering tree, allows sanity check of min_mod_size value

file <- paste(outdir, "clustering_tree_with_modules.pdf", sep = "/")
pdf(file, height = 9, width = 12)
dynamiccolours <- labels2colors(dynamicmods)
table(dynamiccolours)
plotDendroAndColors(genetree, dynamiccolours, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colours")
dev.off()

# Merging modules with similar expression patterns

melist <- moduleEigengenes(datexpr0, colors = dynamiccolours)
mes <- melist$eigengenes
mediss <- 1 - cor(mes)
metree <- hclust(as.dist(mediss), method = "average")
file <- paste(outdir, "clustering_tree_with_merged_modules.pdf", sep = "/")
pdf(file, height = 9, width = 12)
plot(metree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = medissthres, col = "red")
dev.off()
merge <- mergeCloseModules(datexpr0, dynamiccolours, cutHeight = medissthres,
verbose = 3)
mergedcolours <- merge$colors
mergedmes <- merge$newMEs

# Plot a comparison of merged and unmerged modules

file <- paste(outdir, "clustering_tree_compare_modules.pdf", sep = "/")
pdf(file, height = 9, width = 12)
plotDendroAndColors(genetree, cbind(dynamiccolours, mergedcolours),
c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

# Save output for further analyses

modulecolours <- mergedcolours
colourorder <- c("grey", standardColors(50))
modulelabels <- match(modulecolours, colourorder) - 1
mes <- mergedmes
Rfile <- paste(outdir, "modules.RData", sep = "/")
save(mes, modulelabels, modulecolours, dynamiccolours, genetree, tom, file = Rfile)
