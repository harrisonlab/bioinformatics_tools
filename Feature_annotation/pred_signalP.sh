#!/usr/bin/env bash
#SBATCH -J signalP
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4

# set -u
# set -e

CurPath=$PWD
InFile="$1"
if [ -n "$2" ]; then
  SigP_Version="$2"
fi

Source=$(echo $InFile | rev | cut -d "/" -f4 | rev | sed s/_split//g)
Organism=$(echo $InFile | rev | cut -d "/" -f3 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f2 | rev)
InName=$(echo $InFile | rev | cut -d "/" -f1 | rev)
OutFile=$(echo $InName | sed s/.aa//)

if [ -z $SigP_Version ]; then
  WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}/"${Strain}_${InName}_2"
elif [ $SigP_Version == signalp-3.0 ]; then
  WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}/"${Strain}_${InName}_3"
elif [ $SigP_Version == signalp-4.1 ]; then
  WorkDir=$CurPath/${SLURM_JOB_USER}_${SLURM_JOBID}/"${Strain}_${InName}_4"
elif [ $SigP_Version == signalp-5.0 ]; then
  WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}/"${Strain}_${InName}_5"
fi

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile proteins.fa

# Note - Cluster updates have caused SignalP to stop working when run from
# any location other than a /tmp/ drive on a worker node.
# To deal with this, the script now installs signalP in the tempory directory
# associated with each job.
# Oddly, when the program was installed within the SGE variable $TMPDIR it
# didnt run. This seems to be linked to the length of the PATH to the directory
# where signalP is run from.


if [ -z $SigP_Version ]; then
  echo '$2 unset: Running using default behaviour - SignalP-2.0'
  echo "Installing SignalP"
  Tarball=$(ls /home/armita/prog/signalP/signalp-2.0/signalp-2.0.Linux.tar.Z)
  cp $Tarball .
  tar -zxf *.tar.Z
  # sed -i '2s/^/set -x\n/' signalp-2.0/signalp
  sed -i "s&SIGNALP=.*&SIGNALP=$PWD/signalp-2.0&g" signalp-2.0/signalp
  sed -i 's&then AWK=gawk&then AWK=/usr/bin/awk&g' signalp-2.0/signalp
  echo "SignalP installed"
  ls -lh signalp-2.0/syn/
  # echo "Running test"
  # signalp-2.0/signalp -t euk signalp-2.0/test/test.seq > signalp-2.0/tmp.txt
  # echo "test run"
  signalp-2.0/signalp -t euk -f summary -trunc 70 "proteins.fa" > "$OutFile"_sp.tmp
  rm -r signalp-2.0
  # signalp-2.0 -t euk -f summary -trunc 70 "proteins.fa" > "$OutFile"_sp.tmp
  tail -n +5 "$OutFile"_sp.tmp > "$OutFile"_sp.txt
  echo '----------------------------------------------------------------------' >> "$OutFile"_sp.txt
  rm "$OutFile"_sp.tmp
  PathToAnnotateSigP=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  $PathToAnnotateSigP/annotate_signalP2hmm3_v3.pl "$OutFile"_sp.txt "$OutFile"_sp.tab "$OutFile"_sp.aa "$OutFile"_sp_neg.aa "proteins.fa"
  OutDir=$CurPath/gene_pred_vAG/"$Source"_sigP/$Organism/$Strain/split
elif [ $SigP_Version == signalp-3.0 ]; then
  echo "Running using SignalP version: $SigP_Version"
  echo "Installing SignalP"
  Tarball=$(ls /home/armita/prog/signalP/signalp-3.0/signalp-3.0.Linux.tar.Z)
  cp $Tarball .
  tar -zxf *.tar.Z
  sed -i "s&SIGNALP=.*&SIGNALP=$PWD/signalp-3.0&g" signalp-3.0/signalp
  sed -i 's&AWK=nawk&AWK=/usr/bin/awk&g' signalp-3.0/signalp
  sed -i 's&AWK=/usr/bin/gawk&AWK=/usr/bin/awk&g' signalp-3.0/signalp
  echo "SignalP installed"
  ls -lh signalp-3.0/syn/
  echo "Running test"
  signalp-3.0/signalp -t euk signalp-3.0/test/test.seq > signalp-3.0/tmp.txt
  echo "test run"
  signalp-3.0/signalp -t euk -f short -trunc 70 "proteins.fa" > "$OutFile"_sp.txt
  rm -r signalp-3.0
  # $SigP_Version -t euk -f short -trunc 70 "proteins.fa" > "$OutFile"_sp.txt
  echo '----------------------------------------------------------------------' >> "$OutFile"_sp.txt
  PathToAnnotateSigP=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  $PathToAnnotateSigP/sigP_3.0_parser.py --inp_sigP "$OutFile"_sp.txt --out_tab "$OutFile"_sp.tab --out_fasta "$OutFile"_sp.aa --out_neg "$OutFile"_sp_neg.aa --inp_fasta "proteins.fa"
  OutDir=$CurPath/gene_pred/"${Source}_${SigP_Version}"/$Organism/$Strain/split
elif [ $SigP_Version == signalp-4.1 ]; then
  echo "Running using SignalP version: $SigP_Version"
  /home/agomez/scratch/apps/prog/signalp-4.1/$SigP_Version -t euk -f summary -c 70 "proteins.fa" > "$OutFile"_sp.txt
  echo '----------------------------------------------------------------------' >> "$OutFile"_sp.txt
  PathToAnnotateSigP=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Feature_annotation
  $PathToAnnotateSigP/sigP_4.1_parser.py --inp_sigP "$OutFile"_sp.txt --out_tab "$OutFile"_sp.tab --out_fasta "$OutFile"_sp.aa --out_neg "$OutFile"_sp_neg.aa --inp_fasta "proteins.fa"
  OutDir=$CurPath/gene_pred/"${Source}_${SigP_Version}"/$Organism/$Strain/split
  elif [ $SigP_Version == signalp-5.0 ]; then
  echo "Running using SignalP version: $SigP_Version"
  $SigP_Version -org euk -format short -fasta "proteins.fa" -verbose -prefix "$OutFile"
  echo '----------------------------------------------------------------------' >> "$OutFile"_summary.signalp5
  PathToAnnotateSigP=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Feature_annotation
  $PathToAnnotateSigP/sigP_5.0_parser.py --inp_sigP "$OutFile"_summary.signalp5 --out_tab "$OutFile"_sp.tab --out_fasta "$OutFile"_sp.aa --out_neg "$OutFile"_sp_neg.aa --inp_fasta "proteins.fa"
  OutDir=$CurPath/gene_pred/"${Source}_${SigP_Version}"/$Organism/$Strain/split
else
  echo "SigP program not recognised"
fi
rm proteins.fa
mkdir -p $OutDir
mv * $OutDir/.
rm -r $WorkDir

exit
