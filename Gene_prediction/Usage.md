
# Braker

## Requirements

Conda installation

```bash

#Download GeneMark-ES/ET/EP ver 4.57_lic from http://exon.gatech.edu/GeneMark/license_download.cgi. to /home/gomeza/prog/genemark
#Unpack

cd #change to home directory
cp /home/gomeza/prog/genemark/gm_key_64   

conda install numpy
conda install -c anaconda boost
conda install -c bioconda bamtools
conda install -c bioconda braker
conda install -c bioconda augustus
```




```
NOTE: Authors recommend manual installation for the usage of the latets releases. This section will include info of how to manually install braker and augustus.

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





```bash
for Assembly in $(ls N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f6 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/braker2/$Organism/$Strain
AcceptedHits=alignment/star/N.ditissima/R0905/star_aligmentAligned.sortedByCoord.out.bam
#AcceptedHits=alignment/concatenated.bam
GeneModelName="$Organism"_"$Strain"_braker
#rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```
cp /projects/armita/prog/genemark/2018/gm_key_64  ~/.gm_key

--cores 8 \
#--BAMTOOLS_PATH=/home/armita/prog/bamtools/bamtools/bin 

home/armita/prog/genemark/2019/gm_et_linux_64/gmes_petap 

/home/gomeza/miniconda3/envs/gene_pred/bin/braker.pl \
  --GENEMARK_PATH=/home/gomeza/prog/genemark/gmes_linux_64 \
  --BAMTOOLS_PATH=/home/gomeza/miniconda3/envs/gene_pred/bin \
  --overwrite \
  --fungus \
  --gff3 \
  --softmasking on \
  --species=N.ditissima_R0905_v2 \
  --genome=N.ditissima/R0905/polished/repeat_masked/filtered_contigs/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa \
  --bam=alignment/star/N.ditissima/R0905/star_aligmentAligned.sortedByCoord.out.bam 

  braker.pl \
  --GENEMARK_PATH=/home/gomeza/prog/genemark/gmes_linux_64 \
  --AUGUSTUS_CONFIG_PATH=/home/gomeza/prog/Augustus/config \
  --BAMTOOLS_PATH=/home/gomeza/miniconda3/envs/gene_pred/bin \
  --overwrite \
  --fungus \
  --gff3 \
  --softmasking on \
  --species=N.ditissima_R0905_v4 \
  --genome=R0905_contigs_unmasked.fa \
  --bam=concatenated.bam

  ```bash
for Assembly in $(ls R0905_contigs_unmasked.fa); do
OutDir=gene_pred/braker
AcceptedHits=concatenated.bam
GeneModelName=Nd_R0905_braker
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/braker_fungi_v2.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```
cambie el fiterGenesIN_mRNAname.pl
  ```bash
for Assembly in $(ls R0905_contigs_unmasked.fa); do
OutDir=gene_pred/braker_v2
AcceptedHits=concatenated.bam
GeneModelName=Nd_R09052_braker
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/braker_fungi_v2.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

perl /home/gomeza/miniconda3/envs/gene_pred/bin/filterGenesIn_mRNAname.pl good_genes.lst train.gb > train.f.gb 2>errors/filterGenesIn_mRNAname.stderr
