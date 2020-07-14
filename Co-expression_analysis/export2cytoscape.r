#!/usr/bin/Rscript


# Load libaries

library("WGCNA")
library("optparse")

# Import option recommended in WGCNA documentation

options(stringsAsFactors = FALSE)

# Parse arguments

opt_list <- list(
    make_option("--out_dir", type = "character",
    help = "Directory for output to be written to"),
    make_option("--module", type = "character",
    help = "module to export for visualisation in Cytoscape")
    )

opt <- parse_args(OptionParser(option_list = opt_list))
outdir <- opt$out_dir
module <- opt$module

lfile <- paste(outdir, "Cleaned_data.RData", sep = "/")
lnames <- load(file = lfile)
lfile2 <- paste(outdir, "modules.RData", sep = "/")
lnames2 <- load(file = lfile2)

# Set variables for writing out files for Cytoscape

transcripts <- names(datexpr0)
inmodule <- is.finite(match(modulecolours, module))
modtranscripts <- transcripts[inmodule]
modtom <- tom[inmodule, inmodule]
dimnames(modtom) <- list(modtranscripts, modtranscripts)

# Write out files for Cytoscape

edgename_start <- paste("cyt_inp_edges", module, sep = "_")
edgename <- paste(edgename_start, "txt", sep = ".")
nodename_start <- paste("cyt_inp_nodes", module, sep = "_")
nodename <- paste(nodename_start, "txt", sep = ".")

cyt <- exportNetworkToCytoscape(modtom,
edgeFile = paste(outdir, edgename, sep = "/"),
nodeFile = paste(outdir, nodename, sep = "/"),
weighted = TRUE, threshold = 0.02, nodeNames = modtranscripts,
altNodeNames = modtranscripts, nodeAttr = modulecolours[inmodule])
