# Gene feature annotation tools

## Interproscan

Interproscan was used to give gene models functional annotations.

Note: This is a long-running script. As such, these commands were run using 'screen' to allow jobs to be submitted and monitored in the background. This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interproscan jobs was redirected to a temporary output file named interproscan_submission.log.


```bash
	screen -a
	ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
	for Genes in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
  	echo $Genes
  	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
Following interproscan annotation split files were combined using the following commands:

  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  for Proteins in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
  	Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
  	Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
  	echo "$Organism - $Strain"
  	echo $Strain
  	InterProRaw=gene_pred_vAG/interproscan/Ref_Genomes/$Organism/$Strain/raw
  	$ProgDir/append_interpro.sh $Proteins $InterProRaw
  done

```