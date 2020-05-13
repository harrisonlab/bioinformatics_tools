# RepeatMasker

RepeatMasker screens DNA sequences for interspersed repeats and low complexity sequences.

### Requirements

```bash
conda activate general_tools

conda install repeatmodeler
#conda install repeatmasker
#conda install rmblast
conda install -c bioconda transposonpsi

# After installing RepeatMasker, run configure. Bioconda does not do this.
cd /home/USER_ID/miniconda3/envs/general_tools/share/RepeatMasker/

./configure

# Set execution path of tfr, e.g. /home/USER_ID/miniconda3/envs/USER_ENV/bin/trf
# Add search engine. Option 2 - RMBlast will be used
# Set path where rmblastn and makeblastdb are found, e.g. /home/USER_ID/miniconda3/envs/USER_ENV/bin
# 5. Done to exit

#No nucleotide repeat library is included in RepeatMasker.If needed this can be done with the following commands
#cp /home/gomeza/miniconda3/envs/general_tools/share/RepeatMasker/Libraries/RepeatMasker.lib ${CONDA_PREFIX}/share/RepeatMasker/Libraries
#makeblastdb -dbtype nucl -in ${CONDA_PREFIX}/share/RepeatMasker/Libraries/RepeatMasker.lib
```

## Rename contigs

Rename the sequences in assembly fasta file to have simple names.

```bash
ProgDir=/home/gomeza/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.txt
for Assembly in $(ls assembly/miniasm/F.venenatum/WT_minion/racon_10/pilon_10.fasta); do
  OutDir=$(dirname $Assembly)
  $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/WT_miniasm_pilon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
rm tmp.txt
```


## RepeatMasker and transposonPSI


```bash
  ProgDir=/home/gomeza/git_repos/tools/seq_tools/repeat_masking
  BestAssembly=assembly/miniasm/F.venenatum/WT_minion/racon_10/pilon/WT_miniasm_pilon10_renamed.fasta
  OutDir=repeat_masked5/F.venenatum/WT_minion/miniasm
  #sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
  sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir
```