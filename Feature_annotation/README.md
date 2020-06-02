# Gene feature annotation tools

## Interproscan

Interproscan was used to give gene models functional annotations.

Note: This is a long-running script. As such, these commands were run using 'screen' to allow jobs to be submitted and monitored in the background. This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interproscan jobs was redirected to a temporary output file named interproscan_submission.log.


```bash
screen -a
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
for Genes in $(ls gene_pred/N.ditissima/R0905_test/final/final_genes_appended_renamed.pep.fasta); do
echo $Genes
$ProgDir/interproscan.sh $Genes
done 2>&1 | tee -a interproscan_submisison.log








Following interproscan annotation split files were combined using the following commands:

for file in $(ls $SplitDir/*_split*); do
	#Jobs=$(qstat | grep 'run_interp' | grep 'qw' | wc -l)
	#while [ $Jobs -gt 1 ]; do
		#sleep 10
		#printf "."
		#Jobs=$(qstat | grep 'run_interp' | grep 'qw' | wc -l)
	#done
	#printf "\n"
	#echo $file
	#qsub $InterproDir/run_interproscan.sh $file
	sbatch --array=1-30%8 $InterproDir/run_interproscan.sh $file
done


  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  for Proteins in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
  	Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
  	Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
  	echo "$Organism - $Strain"
  	echo $Strain
  	InterProRaw=gene_pred_vAG/interproscan/Ref_Genomes/$Organism/$Strain/raw
  	$ProgDir/append_interpro.sh $Proteins $InterProRaw
  done

for file in $(ls interproscan/*_split*); do
sbatch /home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation/run_interproscan.sh $file
done

/data/scratch/gomeza/prog/Interproscan/interproscan-5.44-79.0/interproscan.sh
nanopolish_makerange.py $Assembly | parallel --results nanopolish.results -P 8 \
nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r $Reads \
--ploidy $Ploidy \
--max-haplotypes 100000 \
--fix-homopolymers \
--min-candidate-frequency 0.2 \
-b $Aligned \
-g $Assembly -t 4
/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation/splitfile_500.py --inp_fasta $InFile | parallel --out_dir $OutDir \
interproscan.sh -goterms -iprlookup -pa -i $IN_NAME


nanopolish_makerange.py $Assembly | parallel --results nanopolish.results -P 8 \
nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r $Reads \
--ploidy $Ploidy \
--max-haplotypes 100000 \
--fix-homopolymers \
--min-candidate-frequency 0.2 \
-b $Aligned \
-g $Assembly -t 4

InterproDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
SplitfileDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation

SplitDir=gene_pred/interproscan/$Organism/$Strain
mkdir -p $SplitDir
InName=$(basename $InFile)
$SplitfileDir/splitfile_500.py --inp_fasta final/final_genes_appended_renamed.pep.fasta | parallel --out_dir interproscan --out_base Inter_split \
/data/scratch/gomeza/prog/Interproscan/interproscan-5.44-79.0/interproscan.sh -goterms -iprlookup -pa -i interproscan/Inter_split_{1}.fa | parallel
```