# Feature analysis

## Blast pipe searches

### Requirements

```bash
conda activate perly_env
```

```bash
for Assembly in $(ls path/to/assembly/or/genes/file/*.fasta); do # Use files with nucleotides
  Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
  echo "$Organism - $Strain"
  Query=path/to/query/fasta
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  sbatch $ProgDir/blast_pipe.sh $Query protein $Assembly
done
```

## Orthofinder


### Requirements

All programs are already installed in the /scratch directory

```
# Add these lines to profile
PATH=${PATH}:/scratch/software/orthomclSoftware-v2.0.9/bin
PATH=${PATH}:/scratch/software/diamond-0.9.27
PATH=${PATH}:/scratch/software/fastme-2.1.5/bin/bin
PATH=${PATH}:/scratch/software/orthomclSoftware-v2.0.9/bin
```

### Typical run

```bash
IsolateAbrv=Fv_vs_Fg_JGI
WorkDir=analysis/orthology/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted
#mkdir -p $WorkDir/goodProteins
#mkdir -p $WorkDir/badProteins  
```

Format fasta files

```bash
Taxon_code=STR1 # Simple and unique taxon identifier
Fasta_file=$(ls path/to/protein/file/final_genes_appended_renamed.pep.fasta)
Id_field=1 # Indicate what field contains the protein ID
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=ST2 # Simple and unique taxon identifier
Fasta_file=$(ls path/to/protein/file/final_genes_appended_renamed.pep.fasta)
Id_field=1 # Indicate what field contains the protein ID
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

.....
```

Using orthofinder

```bash
for IN_DIR in in $(ls -d $WorkDir/formatted) ; do
sbatch orthofinder.sh $IN_DIR $IsolateAbrv
done
```
