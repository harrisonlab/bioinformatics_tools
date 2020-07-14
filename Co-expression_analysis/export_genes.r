#!/usr/bin/Rscript


# Load libraries

library("WGCNA")
library("optparse")

# Import option recommended in WGCNA documentation

options(stringsAsFactors = FALSE)

# Parse arguments

opt_list <- list(
    make_option("--out_dir", type = "character",
    help = "Directory for output to be written to"),
    make_option("--unmerge", type = "character",
    help = "Y or N, Y exports unmerged modules, N does not")
    )

opt <- parse_args(OptionParser(option_list = opt_list))
outdir <- opt$out_dir
unmerge <- opt$unmerge

if (unmerge == "Y"){
} else if (unmerge == "N"){
} else {
    print("Error: Only Y or N are accepted values for the unmerge option,
    please specify.")
    quit(save = "no", status = 1, runLast = FALSE)
}

lfile <- paste(outdir, "Cleaned_data.RData", sep = "/")
lnames <- load(file = lfile)
lfile2 <- paste(outdir, "modules.RData", sep = "/")
lnames2 <- load(file = lfile2)

# Load list of transcript IDs and write out merged modules

transcripts <- names(datexpr0)
for (module in unique(modulecolours)){
    modgenes <- (modulecolours == module)
    modtranscripts <- transcripts[modgenes]
    filename <- paste("Genes_in_", module, ".txt", sep = "")
    mergedir <- paste(outdir, "merged_modules", sep = "/")
    file <- paste(mergedir, filename, sep = "/")
    dir.create(mergedir)
    write.table(as.data.frame(modtranscripts), file = file, row.names = FALSE,
    col.names = FALSE)
}

# Write out unmerged modules

if (unmerge == "Y"){
    for (module in unique(dynamiccolours)){
        modgenes <- (dynamiccolours == module)
        modtranscripts <- transcripts[modgenes]
        filename <- paste("Genes_in_", module, ".txt", sep = "")
        unmergedir <- paste(outdir, "unmerged_modules", sep = "/")
        file <- paste(unmergedir, filename, sep = "/")
        dir.create(unmergedir)
        write.table(as.data.frame(modtranscripts), file = file,
        row.names = FALSE, col.names = FALSE)
    }
} else if (unmerge == "N"){
    print("Unmerged modules not saved.")
}
