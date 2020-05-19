# Gene prediction tools

Tools used in the prediction and annotation of genes

1. BUSCO: Indentification of single-copy orthologs that should be highly conserved among the closely related species. This tool can be used to evaluated genome completeness and for phylogenetic analysis.

2. Braker: A pipeline for automated prediction of protein coding genes using GeneMark-ES/ET and AUGUSTUS. This allow ab initio and gene model training predictions.


## Busco

Indentification of Benchmarking Universal Single-Copy Orthologs in the genome usgin BUSCO databases.

### Requirements

```bash
conda activate olc_assemblers # e.g.

conda install busco

#Busco databases odb9 are in /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9 (recommended)
#New databases, odb10 will be downloaded in /projects/dbBusco/
```

### Typical run

```bash
  for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/WT_miniasm_racon10_renamed.fasta); do
    Strain=WT_minion
    Organism=F.venenatum
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd9
    sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
  done
```


## Braker

Braker is a combination of GeneMark-ET and Augustus used for gene annotation.

### Requirements

Braker requires the instalation of multiple perl libraries, Augustus and Genemark. There are conda packages for Augustus and Braker. However, these packages are often lagging behing the releasess of AUGUSTUS and BRAKER. Therefore manual installation might be necessary.

### Conda installation

```bash
Conda activate gene_pred # e.g.

conda install numpy
conda install -c anaconda boost
conda install -c bioconda bamtools
conda install -c bioconda braker
conda install -c bioconda augustus

# GeneMark-ET will run if a VALID key file resides in your home directory.
cp /home/gomeza/prog/genemark/gm_key_64 ~/.gm_key
```

### Genemark-ET installation (optional)

```
# In case you need/want to install GeneMark in your home directory. If so, braker_fungi.sh script needs to be edited with the new Genemark path.

# Download GeneMark-ES/ET/EP ver 4.57_lic from http://exon.gatech.edu/GeneMark/license_download.cgi. to your home directory, e.g. /home/gomeza/prog/genemark

tar -xvf gmes_linux_64.tar.gz
gunzip gm_key_64.gz

cd #change to home directory
cp /path/to/genemark/gm_key_64 ~/.gm_key  
```

### Manual installation (optional, ongoing...)

```
#NOTE: Authors recommend manual installation for the usage of the latets releases. This section will include info of how to manually install braker and augustus.

#clone braker repo
#git clone https://github.com/Gaius-Augustus/BRAKER.git
#genemark downloaded to /home/gomeza/prog
#tar -xzf gmes_linux_64.tar.gz 
#gunzip gm_key.gz 
#cd (to your home directory)
#mv prog/genemark/gm_key_64 .gm_key
#git clone https://github.com/Gaius-Augustus/Augustus.git

#git clone https://github.com/samtools/htslib.git
#cd htslib
  #autoheader
  #autoconf
  #./configure --prefix=/home/gomeza/Augustus
  #make
    #make install
  #cd ..

  #git clone https://github.com/samtools/bcftools.git
  #cd bcftools
 #autoheader
  #autoconf
  #./configure  --prefix=/home/gomeza/Augustus
  #make
  #make install

#git clone https://github.com/samtools/samtools.git

#cd samtools
#autoheader
  #autoconf -Wno-syntax
  #./configure --prefix=/home/gomeza/Augustus
  #make
#make install
  #cd ..

cpan install YAML
cpan install Hash::Merge
cpan install Logger::Simple
cpan install Parallel::ForkManager
cpan install MCE::Mutex # /home/gomeza/miniconda3/envs/perly_env/bin/perl Makefile.PL -- NOT OK??

cpan install File::HomeDir
#cpan install Scalar::Util::Numeric # CHOCOLATE/Scalar-Util-Numeric-0.40.tar.gz /usr/bin/make -- NOT OK
conda install -c bioconda perl-scalar-util-numeric
conda install cdbtools

cpan install Math::Utils

#perl change_path_in_perl_scripts.pl "/usr/bin/env perl"
#export PERL5LIB=$HOME/prog/BRAKER/scripts:$PERL5LIB
#PATH=${PATH}:$HOME/prog/BRAKER/scripts
```

### Important configuration

Current version of AUGUSTUS and Braker are the same installed using conda packages (braker version 2.1.5 and AUGUSTUS 3.3.3). However, some perl scripts in the conda package of AUGUSTUS are incompatible with this version of braker, producing unexpected output files. These were updated but not included in the conda package yet.

Therefore, replace them in your conda installation before running braker.

```bash
cp /home/gomeza/miniconda3/envs/gene_pred/bin/filterGenesIn_mRNAname.pl /home/USER_ID/miniconda3/USER_ENV/bin
```

### Typical run

```bash
  for Assembly in $(ls path/to/softmasked/assembly/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev) # Edit to set your ouput directory
    Organism=$(echo $Assembly | rev | cut -d '/' -f6 | rev) # Edit to set your ouput directory
    echo "$Organism - $Strain"
    OutDir=gene_pred/braker/$Organism/$Strain
    AcceptedHits=path/to/your/spliced/aligments/files/*_aligmentAligned.sortedByCoord.out.bam # STAR output, see Genome_aligners folder
    #AcceptedHits=alignment/concatenated.bam # Concatenatented alignment files can be used
    GeneModelName="$Organism"_"$Strain"_braker 
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```


## CodingQuarry


### Conda installation

```bash
# A conda environment with python 2.7 is required
# e.g. conda create --name gene_pred_py27 python=2.7

conda install codingquarry


nano ~/.profile
# The environmental variable QUARRY_PATH is set in your profile adding
export QUARRY_PATH="/home/"USER_ID"/miniconda3/envs/"USER_ENV_py27"/opt/codingquarry-2.0/QuarryFiles/QuarryFiles"

. ~/.profile # Refresh your profile
```

### Typical run

```bash
  for Assembly in $(ls path/to/unmasked/genome/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev) # Edit to set your ouput directory
    Organism=$(echo $Assembly| rev | cut -d '/' -f6 | rev) # Edit to set your ouput directory
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain/
    mkdir -p $OutDir
    GTF=path/to/RNAseq/alignment/assembly/*.gtf # GFT file from stringtie/cufflinks output. See Genome-guided_assemblers scripts
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    sbatch $ProgDir/codingquarry.sh $Assembly $GTF $OutDir
  done
```