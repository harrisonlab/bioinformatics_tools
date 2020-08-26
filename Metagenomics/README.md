# Metagenomics analysis tools

## Centrifuge

Classification of DNA sequences from microbial samples. 

### Requirements 

```bash
# Add centrifuge to your PATH environment variable.
nano ~/.profile
# Copy the following path
PATH=${PATH}:/scratch/software/centrifuge
# Save and refresh your profile using
. ~/.profile
```

Build index instructions

```bash
# This is not needed unless you are using a different database
# If you do, log into any long node
screen -a
srun --partition long --mem 20G --cpus-per-task 10 --pty bash
# Download NCBI taxonomy
centrifuge-download -o taxonomy taxonomy

# Download ALL complete archaeal, bacterial and viral genomes from NCBI
centrifuge-download -o library -m -d "archaea,bacteria,viral" refseq > seqid2taxid.map
centrifuge-download -o library -m -d "viral" refseq > seqid2taxid.map # only viruses
# Concatenate all sequences
cat prog/centrifuge/library/viral/*.fna > prog/centrifuge/library/input-viral-sequences.fna
# Build centrifuge index with 4 threads
centrifuge-build -p 4 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp input-sequences.fna allvirus

# You can download only genome of interest from the NCBI website using the different filter options. This was done for plant viruses only.
centrifuge-build -p 4 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp plantvirus_sequences.fasta plantvirus
```

### Run centrifuge


```bash
# Log into an interactive session and working node
screen -a
srun --partition long --mem 20G --cpus-per-task 10 --pty bash

# Run centrifuge
OutDir=path/to/output/dir
centrifuge -p 4 -x /data/scratch/gomeza/prog/centrifuge/plantvirus -t -q -1 path/to/unmapped/mate1 -2 path/to/unmapped/mate2 --phred33 --report-file $OutDir/centrifuge_report.tsv -S $OutDir/centrifuge_results.txt 
# --min-hitlen <integer> option can be used to set a minimum length of partial hits. Default 22.

# Create a Kraken-style report
centrifuge-kreport -x /data/scratch/gomeza/prog/centrifuge/plantvirus centrifuge_results.txt > centrifuge_krakened.txt
# --min-score <integer> option set minimum score for reads to be counted
# --min-length <integer> option set minimum alignment length to the read
```

Results can be visualised using Pavian. https://fbreitwieser.shinyapps.io/pavian/

