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
ProgDir=home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
touch tmp.txt
for Assembly in $(ls path/to/assembly/*.fasta); do
  OutDir=$(dirname $Assembly)
  $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/WT_miniasm_pilon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
rm tmp.txt
```


## RepeatMasker and transposonPSI


```bash
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
  BestAssembly=path/to/assembly/*renamed.fasta
  OutDir=repeat_masked/$Organism/$Strain/$AssemblyMethod
  sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
  sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir
```


The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.


```bash
for File in $(ls repeat_masked/F.venenatum/WT_minion/miniasm/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/F.venenatum/WT_minion/miniasm/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```
