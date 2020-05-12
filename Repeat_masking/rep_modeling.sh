#!/usr/bin/env bash
#SBATCH -J repeatmasker
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=24G
#SBATCH --cpus-per-task=24

# This script uses repeatmodeler and repeatmasker to mask Interspersed repeats
# and low complexity regions within the genome. Firstly, repeatmodeler identifies
# repeat element boundaries and relationships within repeat families. The repeats
# identified within the genome are provided to repeatmasker, which uses this data
# along with it's own repeat libraries to identify these repetitive regions and
# perform masking. Masking is done at 3 levels:
# Hardmasking = repetitive sequence is replaced with N's.
# Softmasking = repetitive sequence is converted to lower case.
# Ignoring low-complexity regions = only interspersed repetitive elements are masked.

# Edition notes
# Note 1 - The latest version of RepeatModeler (RepeatModeler 2.0.1) copies the final classified results to the families.stk and consensus.fa.
# RepeatClassifier is needed in this version to produced consensi.fa.classified files.
# Note 2 - No nucleotide repeat library is included in RepeatMasker. It is recommended to download a database from RebBase but a licence is needed.
# The RepBase v19.09 was copied to /home/gomeza/RepMaskerDB

Usage="sbatch rep_modeling.sh <assembled_contigs.fa> [<output_directory>]"

InFile=$1
Organism=$(echo $InFile | rev | cut -d "/" -f4 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f3 | rev)
Assembly=$(echo $InFile | rev | cut -d "/" -f2 | rev)

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

if [ $2 ]; then
  OutDir=$CurPath/$2
else
  OutDir=$CurPath/repeat_masked/$Organism/$Strain/"$Assembly"_repmask
fi

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile "$Strain"_contigs_unmasked.fa
BuildDatabase -name "$Strain"_RepMod "$Strain"_contigs_unmasked.fa
RepeatModeler -pa 16 -database "$Strain"_RepMod

RepeatClassifier -consensi RM_*.*/consensi.fa -stockholm RM_*.*/families.stk

# hardmask
RepeatMasker -gff -pa 16 -lib RM_*.*/consensi.fa.classified "$Strain"_contigs_unmasked.fa
mv "$Strain"_contigs_unmasked.fa.cat.gz "$Strain"_contigs_hardmasked.fa.cat.gz
mv "$Strain"_contigs_unmasked.fa.masked "$Strain"_contigs_hardmasked.fa
mv "$Strain"_contigs_unmasked.fa.out "$Strain"_contigs_hardmasked.out
mv "$Strain"_contigs_unmasked.fa.out.gff "$Strain"_contigs_hardmasked.fa.out.gff
mv "$Strain"_contigs_unmasked.fa.tbl "$Strain"_contigs_hardmasked.tbl
grep -v '#' "$Strain"_contigs_hardmasked.fa.out.gff > "$Strain"_contigs_hardmasked.gff

# softmask
RepeatMasker -xsmall -gff -pa 16 -lib RM_*.*/consensi.fa.classified "$Strain"_contigs_unmasked.fa
mv "$Strain"_contigs_unmasked.fa.cat.gz "$Strain"_contigs_softmasked.fa.cat.gz
mv "$Strain"_contigs_unmasked.fa.masked "$Strain"_contigs_softmasked.fa
mv "$Strain"_contigs_unmasked.fa.out "$Strain"_contigs_softmasked.out
mv "$Strain"_contigs_unmasked.fa.out.gff "$Strain"_contigs_softmasked.fa.out.gff
mv "$Strain"_contigs_unmasked.fa.tbl "$Strain"_contigs_softmasked.tbl
grep -v '#' "$Strain"_contigs_softmasked.fa.out.gff > "$Strain"_contigs_softmasked.gff

# don't mask low-complexity or simple-repeat sequences, just transposons
RepeatMasker -nolow -gff -pa 16 -lib RM_*.*/consensi.fa.classified "$Strain"_contigs_unmasked.fa
mv "$Strain"_contigs_unmasked.fa.cat.gz "$Strain"_contigs_transposonmasked.fa.cat.gz
mv "$Strain"_contigs_unmasked.fa.masked "$Strain"_contigs_transposonmasked.fa
mv "$Strain"_contigs_unmasked.fa.out "$Strain"_contigs_transposonmasked.out
mv "$Strain"_contigs_unmasked.fa.out.gff "$Strain"_contigs_transposonmasked.fa.out.gff
mv "$Strain"_contigs_unmasked.fa.tbl "$Strain"_contigs_transposonmasked.tbl
grep -v '#' "$Strain"_contigs_transposonmasked.fa.out.gff > "$Strain"_contigs_transposonmasked.gff

mkdir -p $OutDir
rm -r $WorkDir/RM*
cp -r $WorkDir/* $OutDir/.

rm -r $WorkDir
