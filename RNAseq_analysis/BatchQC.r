# BatchQC 

```r
library(BatchQC)
library(readr) #If tsv needed

# Metadata file with sample IDs in the first column and suspected batch variables
metadf <- read.csv("Vd_clocks_RNAseq_design.csv")

# Gene by sample matrix with gene IDs in the first column and sample IDs as column headers
exprdata <- read.csv("TPM_quantile.csv")


batch = metadf$Experiment

condition = metadf$Condition

 
batchQC(dat=exprdata, batch=batch, condition=condition, 
report_file="batchqc_report.html", report_dir=".",
report_option_binary="111111111",
view_report=TRUE, interactive=TRUE, batchqc_output=TRUE)
```