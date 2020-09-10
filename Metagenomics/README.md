# Metagenomics analysis tools

## Kraken

### Requirements 

```bash
# Add centrifuge to your PATH environment variable.
nano ~/.profile
# Copy the following path
PATH=${PATH}:/scratch/software/kraken2/kraken2-master
# Save and refresh your profile using
. ~/.profile
```
 
```bash
# This is not needed unless you are using a different database
# If you do, log into any long node
screen -a
srun --partition long --mem 20G --cpus-per-task 10 --pty bash
# Download NCBI taxonomy
kraken2-build --download-taxonomy --db plantvirusesDB/
# Add sequences to library. See manual for the installation of other databases.
kraken2-build --add-to-library plantviruses.fasta --db plantvirusesDB
# Build database
kraken2-build --build --db plantvirusesDB/
```

### Run centrifuge


```bash
for Read1 in $(ls alignment/star/*/star_aligmentUnmapped.out.mate1); do
    Read2=$(echo $Read1 | sed 's/mate1/mate2/g')
    echo $Read1
    echo $Read2
    Sample=$(echo $Read1 | rev | cut -f2 -d '/' | rev)
    Database=plantvirus
    OutDir=analysis/Centrifuge/$Sample
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
    sbatch $ProgDir/centrifuge.sh $Database $Read1 $Read2 $OutDir
done
```

```bash
# Log into an interactive session and working node
screen -a
srun --partition long --mem 20G --cpus-per-task 10 --pty bash

# Run kraken2
OutDir=path/to/output/dir
kraken2 --db /data/scratch/gomeza/prog/kraken2/plantvirusesDB --threads 4 --paired --classified-out cseqs#.fq path/to/unmapped/mate1 path/to/unmapped/mate2 --output $OutDir/KrakenResults.txt
```

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

##### Build index instructions

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

centrifuge-download -o library -m -t "NC_024710	" refseq > seqid2taxid.map
dustmasker -in plantvirus_sequences.fasta -infmt fasta -out sequences -outfmt fasta
centrifuge-build -p 4 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp plantvirus_masked.fasta plantvirus

```
```bash
ProgDir=/home/gomeza/git_repos/tools/seq_tools/repeat_masking
BestAssembly=plantvirus_sequences.fasta
OutDir=repeat_masked/
sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir
```

```bash
for File in $(ls repeat_masked/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```
centrifuge-build -p 4 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp repeat_masked/plantvirus_sequences.fasta_contigs_hardmasked.fa plantvirushard

### Run centrifuge


```bash
for Read1 in $(ls alignment/star/*/star_aligmentUnmapped.out.mate1); do
    Read2=$(echo $Read1 | sed 's/mate1/mate2/g')
    echo $Read1
    echo $Read2
    Sample=$(echo $Read1 | rev | cut -f2 -d '/' | rev)
    Database=plantvirus
    OutDir=analysis/Centrifuge/$Sample
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
    sbatch $ProgDir/centrifuge.sh $Read1 $Read2 $Database $OutDir
done
```

Results can be visualised using Pavian. https://fbreitwieser.shinyapps.io/pavian/



### Run centrifuge


```bash
for Read1 in $(ls Charlotte/star_aligmentUnmapped.out.mate1); do
Read2=$(echo $Read1 | sed 's/mate1/mate2/g')
echo $Read1
echo $Read2
Sample=$(echo $Read1 | rev | cut -f2 -d '/' | rev)
Database=plantvirusesDB
OutDir=analysis/Kraken/$Sample
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
sbatch $ProgDir/kraken.sh $Read1 $Read2 $Database $OutDir
done
```