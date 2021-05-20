# Commands for phylogenetic analyses 

## BUSCO genes

For BUSCO identification see Gene_prediction folder

### Requirements

```bash
conda activate general_tools
# Gene aligner. Version=v7.475
conda install -c bioconda mafft
# Trim poor aligned sequences. Version=1.4.1
conda install -c bioconda trimal
# Randomized Axelerated Maximum Likelihood. Version=8.2.12
conda install -c bioconda raxml
# Process phylogenetic trees. Version=1.6.4
conda install newick_utils
```

## Find single copy busco genes

Create a list of all BUSCO IDs

```bash
    OutDir=analysis/popgen/busco_phylogeny
    mkdir -p $OutDir
    BuscoDb="sordariomycetes_odb10"
    ls -1 /projects/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt
```

```bash
# For BUSCO version 4
screen -a
srun --partition long --mem-per-cpu 10G --cpus-per-task 10 --pty bash
# Create a folder for each busco gene
mkdir temp_busco
printf "" > analysis/popgen/busco_phylogeny/single_hits.txt
  for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
  echo $Busco
  OutDir=analysis_VP/popgen/busco_phylogeny_2/$Busco
  mkdir -p $OutDir
  # Move all single copy genes to each folder and rename gene headers
    for Fasta in $(ls gene_pred/busco/$Organism/$Strain/*/*/*/*/single_copy_busco_sequences/$Busco*.fna); do
      Strain=$(echo $Fasta | rev | cut -f7 -d '/' | rev)
      Organism=$(echo $Fasta | rev | cut -f8 -d '/' | rev)
      FileName=$(basename $Fasta)
      contig=$(cat $Fasta | grep '>' | sed 's/ <unknown description>//g' | sed 's/>//g')
      echo ">$Busco:$Strain:$contig" > temp_busco/"$Busco"_"$Strain"_new_names.txt
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools
      python $ProgDir/replace_fasta_records.py -i $Fasta -r temp_busco/"$Busco"_"$Strain"_new_names.txt -o $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
      rm temp_busco/"$Busco"_"$Strain"_new_names.txt
    #cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
    done
  # Create fasta file containing all busco for alignment
  cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
  SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
  printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny/single_hits.txt
  done
rm -r temp_busco
```

```bash
# For BUSCO version 3 only (deprecated)
# log in a worker node
screen -a
srun --partition long --mem-per-cpu 10G --cpus-per-task 10 --pty bash
# Create a folder for each busco gene
printf "" > analysis/popgen/busco_phylogeny/single_hits.txt
  for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
  echo $Busco
  OutDir=analysis/popgen/busco_phylogeny/$Busco
  mkdir -p $OutDir
# Move all single copy genes to each folder
  for Fasta in $(ls gene_pred/busco/$Organism/$Strain/*/*/*/*/single_copy_busco_sequences/$Busco*.fna); do
  Strain=$(echo $Fasta | rev | cut -f5 -d '/' | rev)
  Organism=$(echo $Fasta | rev | cut -f6 -d '/' | rev)
  FileName=$(basename $Fasta)
  cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
  done
# Create fasta file containing all busco for alignment
  cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
  SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
  printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny/single_hits.txt
done
```

```bash
# Check for multiple hits and remove
less analysis/popgen/busco_phylogeny/single_hits.txt | sort -k2 -n
# e.g.
rm gene_pred/busco/N.ditissima/*/*/*/single_copy_busco_sequences/EOG093318S0*
rm gene_pred/busco/N.ditissima/*/*/*/single_copy_busco_sequences/EOG093305B4*
rm gene_pred/busco/N.ditissima/*/*/*/single_copy_busco_sequences/EOG0933010Y*
```

```bash
# If all isolates have a single copy of a busco gene, move the appended fasta to a new folder
OutDir=analysis/popgen/busco_phylogeny/alignments
mkdir -p $OutDir
OrganismNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | cut -f2 | sort -nr | head -n1)
for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
echo $Busco
HitNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | grep "$Busco" | cut -f2)
if [ $HitNum == $OrganismNum ]; then
  cp analysis/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta $OutDir/.
fi
done
```

## Gene alignments

```bash
# Submit alignment for single copy busco genes with a hit in each organism
AlignDir=analysis/popgen/busco_phylogeny/alignments
CurDir=$PWD
cd $AlignDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics
squeue $ProgDir/mafft.sh
cd $CurDir
```

## Nucleotide diversity (optional)

For closely related organisms, identify genes with high nucleotide diversity (Pi) and average number of pairwise differences, medium number of segregating sites (avoid alignments with low homology and lots of phylogenetically uninformative singletons).

For analyses involving cross-species comparisons involving highly diverged sequences with high nucleotide diversity (e.g. [0.1<Pi<0.4]), looking for genes with the lowest number of segregating sites.

```bash
    cd analysis/popgen/busco_phylogeny/alignments
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics
    python $ProgDir/calculate_nucleotide_diversity.py "*aligned.fasta"
```

## Trim poor alignments

Trimming sequence alignments using Trim-Al. Note - automated1 mode is optimised for ML tree reconstruction

```bash
  OutDir=analysis/popgen/busco_phylogeny/trimmed_alignments
  mkdir -p $OutDir
  for Alignment in $(ls analysis/popgen/busco_phylogeny/alignments/*_appended_aligned.fasta); do
    TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
    echo $Alignment
    trimal -in $Alignment -out $OutDir/$TrimmedName -keepheader -automated1
  done
```

## Randomized Axelerated Maximum Likelihood

```bash
# Edit header name keeping BUSCO name and isolate name using sed
cd analysis/popgen/busco_phylogeny/trimmed_alignments
# e.g.
sed -i 's/:contig_.*//g' *_appended_aligned_trimmed.fasta
sed -i 's/:LD.*//g' *_appended_aligned_trimmed.fasta
sed -i 's/:NODE.*//g' *_appended_aligned_trimmed.fasta
```

```bash
screen -a
# sleep option will submit a job every 10s to the short queue
    for Alignment in $(ls analysis/popgen/busco_phylogeny/trimmed_alignments/*aligned_trimmed.fasta); do
        sleep 10s
        Prefix=$(basename $Alignment | cut -f1 -d '_')
        OutDir=analysis/popgen/busco_phylogeny/RAxML/$Prefix
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics
        squeue $ProgDir/RAxML.sh $Alignment $Prefix $OutDir
    done
```

## Astral

Run Astral to build a consensus phylogeny from a collective set of "best phylogenies" from each BUSCO locus.
Note - "Recent versions of ASTRAL output a branch support value even without bootstrapping. Our analyses have revealed that this form of support is more reliable than bootstrapping (under the conditions we explored). Nevertheless, you may want to run bootstrapping as well."
Tutorial tips: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial-template.md

```bash
screen -a 
# Edit if necessary. This will run on compute10
srun --partition long --mem-per-cpu 10G --cpus-per-task 24 --pty bash 

OutDir=analysis/popgen/busco_phylogeny/ASTRAL
mkdir -p $OutDir

#Concatenate best trees
Name=Name4yourPhylogeny
#cat analysis_VP/popgen/busco_phylogeny/RAxML/*/RAxML_bestTree.*  | sed -r "s/CTG.\w+:/:/g" > $OutDir/Nd_phylogeny.appended2.tre
cat analysis/popgen/busco_phylogeny/RAxML/*/RAxML_bestTree.*  > $OutDir/"$Name"_phylogeny.appended.tre

# Contract low support brances (below 10% bootsrap support)
nw_ed $OutDir/"$Name"_phylogeny.appended.tre 'i & b<=10' o > $OutDir/"$Name"_phylogeny.appended.trimmed.tre

# Calculate combined tree
ProgDir=/scratch/software/ASTRAL/ASTRAL-5.7.1/Astral
java -jar $ProgDir/astral.5.7.1.jar -i $OutDir/"$Name"_phylogeny.appended.tre -o $OutDir/"$Name"_phylogeny.consensus.tre | tee 2> $OutDir/"$Name"_phylogeny.consensus.log
# Score the resulting tree
java -jar $ProgDir/astral.5.7.1.jar -q $OutDir/"$Name"_phylogeny.consensus.tre -i $OutDir/"$Name"_phylogeny.appended.tre -o $OutDir/"$Name"_phylogeny.consensus.scored.tre 2> $OutDir/"$Name"_phylogeny.consensus.scored.log
```

Manual edition of the final consensus tree is needed 

```
Step 1: Download consensus tree to local machine

Step 2: Import into geneious and export again in newick format to get around polytomy branches having no branch length.

Step 3: Terminal branch lengths are meanlingless from ASTRAL and should all be set to an arbitrary value. This will be done by geneious (set to 1), but it also introduces a branch length of 2 for one isolate that needs to be corrected with sed
```
```bash
cat Alt_phylogeny.consensus.scored.geneious.tre | sed 's/:2/:1/g' > Alt_phylogeny.consensus.scored.geneious2.tre
```

## Plot best scored tree

GGtree was used to make a plot. Tutorial tips: https://bioconnector.org/r-ggtree.html

R version > 4.0

```r
setwd("/data/scratch/gomeza/")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)
library(ggtree)
library(phangorn)
library(treeio)





tree <- read.tree("/Users/armita/Downloads/Aalt/ASTRAL/expanded/Alt_phylogeny.consensus.scored.geneious2.tre")

t<-ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()

mydata <- read.csv("/Users/armita/Downloads/Aalt/ASTRAL/traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$Isolate
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]


# Format nodes by values
nodes <- data.frame(t$data)
#nodes <- nodes[!nodes$isTip,]
nodes$label <- as.numeric(nodes[nodes$label,])
as.numeric(nodes$label)
#nodes$label[nodes$label < 0.80] <- ''
nodes$support[nodes$isTip] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label > 0.80)] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label < 0.80)] <- 'unsupported'
nodes$support[(!nodes$isTip) & (nodes$label == '')] <- 'supported'
t <- t + aes(linetype=nodes$support)
nodes$label[nodes$label > 0.80] <- ''
t <- t + geom_nodelab(data=nodes, size=2, hjust=-0.05) # colours as defined by col2rgb



#51 is the node of your outgroup?
tree$edge.length[tree$edge.length == 1] <- 0
tree$edge.length[51] <- 0


t <- ggtree(tree, aes(linetype=nodes$support)) # Core tree


# Adjust terminal branch lengths:
branches <- t$data

branches <- t$data
tree$edge.length[branches$isTip] <- 1.0
# tree$edge.length[tree$edge.length == 1] <- 0
# t <- ggtree(tree, aes(linetype=nodes$support))
#Tree <- branches$branch.length

t <- t + geom_treescale(offset=-1.0, fontsize = 3) # Add scalebar
# t <- t + xlim(0, 0.025) # Add more space for labels


# Colouring labels by values in another df
t <- t %<+% mydata # Allow colouring of nodes by another df
#t <- t + geom_tiplab(aes(color=Source), size=3, hjust=0) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

tips <- data.frame(t$data)
tips$label <- tips$ID
t <- t + geom_tiplab(data=tips, aes(color=Source), size=3, hjust=0, align=T, offset = +0.1) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

# Add in a further set of labels
tree_mod <- data.frame(t$data)
tree_mod$label <- tips$pathotype
t <- t + geom_tiplab(data=tree_mod, aes(label=label, color=Source), align=T, linetype = NULL, size=3, hjust=0, offset = +5.0) +
scale_color_manual(values=c("gray39","black"))

tips$MAT <- factor(tips$MAT)
# t <- t + geom_tippoint(data=tips, aes(shape=MAT), size=2)
t <- t + geom_tiplab(data=tips, aes(label=MAT, color=Source), align=T, linetype = NULL, size=3, hjust=0, offset = +3.5) +
scale_color_manual(values=c("gray39","black"))



# Annotate a clade with a bar line
# t <- t + geom_cladelabel(node=42, label='sect. Alternaria', align=T, colour='black', offset=-1.5)
# t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=-4.5)
# t <- t + geom_cladelabel(node=51, label='tenuissima clade', align=T, colour='black', offset=-4.5)
# t <- t + geom_cladelabel(node=45, label='arborescens clade', align=T, colour='black', offset=-4.5)
t <- t + geom_cladelabel(node=43, label='sect. Alternaria', align=T, colour='black', offset=9.5)
t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=6.5)
t <- t + geom_cladelabel(node=46, label='tenuissima clade', align=T, colour='black', offset=6.5)
t <- t + geom_cladelabel(node=65, label='arborescens clade', align=T, colour='black', offset=6.5)
t <- t + geom_cladelabel(node=65, label='', colour='NA', offset=17.5)

# Save as PDF and force a 'huge' size plot
# t <- ggsave("expanded/Fig3_busco_phylogeny.pdf", width =30, height = 30, units = "cm", limitsize = FALSE)
t <- ggsave("expanded/Fig3_busco_phylogeny.tiff", width =30, height = 30, units = "cm", limitsize = FALSE)

````








```
Visually inspect the alignments of selected genes (genes_selected_for_phylogeny.txt) to be used in
constructing the phylogenies and trim them as necessary in MEGA7.
Copy the relevant trimmed alignment FASTA files into

```bash
  # mkdir $CurDir/beast_runs/candidates/select/trimmed
```


##PartitionFinder (nucleotide sequence evolution model)

```bash
cd analysis/popgen/busco_phylogeny/phylogeny

config_template=/home/sobczm/bin/PartitionFinder1.1.1/partition_finder.cfg
ct=$(basename "$config_template")

mkdir NEXUS

# prepare directory for PartitionFinder run:
for f in $(ls *fasta); do
sed -i 's/:/_/g' $f
c="$(cat $f | awk 'NR%2==0' | awk '{print length($1)}' | head -1)"
p="${f%.fasta}.phy"
n="${f%.fasta}.NEXUS"
dir="${f%.fasta}"

mkdir $dir
cp $config_template $dir/.

# Substitute the name of the alignment file and the sequence length in the config file to become correct for the current run.
sed -i 's,^\(alignment = \).*,\1'"$p;"',' $dir/$ct
sed -i 's,^\(Gene1_pos1 = \).*,\1'"1-$c\\\3;"',' $dir/$ct
sed -i 's,^\(Gene1_pos2 = \).*,\1'"2-$c\\\3;"',' $dir/$ct
sed -i 's,^\(Gene1_pos3 = \).*,\1'"3-$c\\\3;"',' $dir/$ct

# Convert FASTA to phylip for the Partition Finder run
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
$ProgDir/fasta2phylip.pl $f>$p
mv $p $dir

# Convert FASTA to NEXUS for the BEAST run
$ProgDir/Fasta2Nexus.pl $f>$dir/$n

#Problems running PartitionFinder on the cluster. May have to be run locally on your Mac or Windows machine.
# qsub $ProgDir/sub_partition_finder.sh $dir
done
```

Partition finder wasnt run on the cluster. As such fasta alignment files were
downloaded to the local machine where partitionfinder was run
patritionfinder2 was downloaded from:
http://www.robertlanfear.com/partitionfinder/

and the anaconda libraries to support it were downloaded from:
https://www.continuum.io/downloads#macos


copy the fasta files and the partitionfinder config files to
your local computer

```bash
cd Users/armita/Downloads
scp -r cluster:/home/groups/harrisonlab/project_files/idris/analysis/popgen/busco_phylogeny/phylogeny .
```

Alignments were loaded into Geneious where they were visualised and manually sorted into
three categories:
* Good - All sequences present no trimming needed
* Trim - All sequences present short regions may need trimming from the beginning / end of the alignment before use in phylogenetics
* Bad - a region of one or more sequences is missing or the sequences / alignment is not appropriate for phylogenetics

These alignments were then exported from Geneious into the following folders:

```bash
cd Users/armita/Downloads/phylogeny
mkdir good_alignments
mkdir trim_alignments
mkdir bad_alignments
```

Alignments within the "good alignments" directory were taken forward for further
analysis

```bash
  for Dir in $(ls -d *_alignments); do
    for Alignment in $(ls $Dir/*_appended_aligned.phy); do
      Prefix=$(echo $Alignment | cut -f2 -d '/' | sed 's/.phy//g')
      echo $Prefix
      cp $Prefix/$Prefix.NEXUS $Dir/$Prefix/.
      cp -r $Prefix $Dir/.
      /Users/armita/anaconda2/bin/python ../partitionfinder-2.1.1/PartitionFinder.py $Dir/$Prefix --no-ml-tree --force-restart
    done
  done > log.txt
```


Upload partition models back to the cluster:

```bash
ClusterDir=/home/groups/harrisonlab/project_files/idris/analysis/popgen/busco_phylogeny/phylogeny
scp -r bad_alignments cluster:$ClusterDir/.
```


## Preparing to run BEAST


Using trimmed FASTA alignments and nucleotide substitution models identified with PartitionFinder:
create an XML input file using BEAUTi, with StarBeast template.

Prepare a 30 loci dataset, in addition to a 5 loci subset to compare convergence.

Run after qlogin into a worker node (BEAST does not find BEAGLE libraries when using qsub -
as the BEAST package is quite fiddly, may troubleshoot it later when necessary.

StarBeast settings used here:
* Substitution rate: default HKY
* Strict clock
* Species Tree Population Size: Linear with constant root
* Yule prior on species tree
* Chain length: 300 million (this may vary, change run convergence with Tracer during the run to establish the number of iterations required
* Tracer: /home/sobczm/bin/beast/Tracer_v1.6/bin/tracer
some runs may never converge)
* Store every: 10000

```bash

cd /home/groups/harrisonlab/project_files/idris


for File in $(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/analysis/best_scheme.txt); do
Busco=$(echo $File | cut -f6 -d '/' | cut -f1 -d '_')
Model=$(cat $File | grep -A1 'Best Model' | tail -n1 | cut -f2 -d '|')
printf "$Busco\t$Model\n"
done

# Edit NEXUS files:
for Nexus in $(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/*_appended_aligned.NEXUS); do
  sed -i -r "s/^.*_P\./P./g" $Nexus
  sed -i -r "s/_contig.*\t/\t/g" $Nexus
  sed -i -r "s/_NODE.*\t/\t/g" $Nexus
done

# OUtputs of partitionfinder were used to set models
# of DNA evolution in Beauti, as described on:
# http://www.robertlanfear.com/partitionfinder/faq/#toc-beast
# CHain length was modified from 10000000 to 500000000 as determined
# by a first run of beast where tracer reported the estimated sasmple size to be below 100 (3) - increase by 50 fold.

# Run Beauti
NexusFiles=$(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/*.NEXUS | sed -e 's/^/ -nex /g' | tr -d '\n')
OutFile=$(echo $Nexus | sed 's/.NEXUS/.xml/g')
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/beauti -template StarBeast.xml $NexusFiles




qlogin -pe smp 8
InXML=analysis/popgen/busco_phylogeny/phylogeny/Pcac_beauti_starBEAST2.xml
OutDir=$(dirname $InXML)"/BEAST4"
mkdir -p $OutDir
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/beast -threads 8 -prefix $OutDir $InXML > $OutDir/log.txt
# java -Djava.library.path="C:\Program Files (x86)\Common Files\libhmsbeagle-1.0" -jar "/BEAST175/lib/beast.jar"

#After the run, check convergence with Tracer, summarise the final tree with TreeAnnotator
for Tree in $(ls $OutDir/*.trees); do
BurnIn=10 # percentage of states to be considered as burnin
SumTree=$(echo $Tree | sed 's/.trees/_summary.tree/g')
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/treeannotator -heights median -burnin $BurnIn $Tree $SumTree
done

#Visualise and beautify the final tree (suffix "summary") with FigTree
FigTree=/home/sobczm/bin/FigTree_v1.4.2/bin/figtree
$FigTree
```
