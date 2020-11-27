# Genome synteny plots

## Circos plot

This program is used to convert fasta files into input format for circos

### Hg199 miniasm vs R0905 canu

```bash
  OutDir=analysis/circos/Hg_vs_R9_vAG
  mkdir -p $OutDir
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis

  Hg199_genome=$(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/Hg199_pilon10_renamed.fasta)
  $ProgDir/fasta2circos.py --genome $Hg199_genome --contig_prefix "Hg_" > $OutDir/Hg199_genome.txt

  R0905_genome=$(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/R0905_pilon10_renamed.fasta)
  $ProgDir/fasta2circos.py --genome $R0905_genome --contig_prefix "R9_" > $OutDir/R0905_genome.txt

  cat $OutDir/Hg199_genome.txt > $OutDir/Hg199_R0905_genome.txt
  tac $OutDir/R0905_genome.txt >> $OutDir/Hg199_R0905_genome.txt
```

Identify Telomere repeats:
Telomeric repeats were identified in assemblies

```bash
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/Hg199_pilon10_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/telomere_vAG/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
    $ProgDir/annotate_telomeres.py --fasta $Assembly --out $OutDir/telomere_hits
  done

  for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/R0905_pilon10_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/telomere_vAG/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
    $ProgDir/annotate_telomeres.py --fasta $Assembly --out $OutDir/telomere_hits
  done
```

Telomere locations on contigs:

```bash
  OutDir=analysis/circos/Hg_vs_R9_vAG
  cat analysis/telomere_vAG/N.ditissima/Hg199/telomere_hits_circos.txt | sed 's/contig/Hg_contig/g' | sort -k3 -n -t'_' > $OutDir/Hg_vs_R9_telomere_hits.txt
  cat analysis/telomere_vAG/N.ditissima/R0905/telomere_hits_circos.txt  | sed 's/contig/R9_contig/g' | sort -k3 -n -t'_' >> $OutDir/Hg_vs_R9_telomere_hits.txt
```
```bash
  OutDir=analysis/circos/Hg_vs_R9_vAG
  Coords=$(ls analysis/genome_alignment/mummer/N.ditissima/Hg199/Hg199_vs_R0905_vAG/Hg199_vs_R0905_vAG_coords.tsv)
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
  $ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id Hg --ref_id R9 > $OutDir/Hg_vs_R9_links.txt
  cat $OutDir/Hg_vs_R9_links.txt > $OutDir/Hg_vs_R9_links_edited.txt
```

A file showing contig orientations was made:
```bash
  cat $OutDir/Hg199_R0905_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/Hg_contig_order.txt
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
  $ProgDir/find_contig_orientation.py --links_file $OutDir/Hg_vs_R9_links_edited.txt > $OutDir/Hg_vs_R9_contig_orientation.txt
```

Contig order was selected by taking the first line of that file and then also taking the reversed order of contigs using the command:

```bash
  cat $OutDir/Hg_vs_R9_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
  cat $OutDir/Hg_vs_R9_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1
  cat $OutDir/Hg199_R0905_genome.txt | grep 'Hg' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Hg/, Hg/g'
  cat $OutDir/Hg199_R0905_genome.txt | grep 'R9' | cut -f3 -d ' ' | tr -d '\n' | sed 's/R9/, R9/g' >> tmp.txt
  echo "Order of unseen Hg contigs and remaining R9 contigs"
  cat $OutDir/Hg199_R0905_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/R9/, R9/g' | sed 's/Hg/, Hg/g'
```

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos/Hg_vs_R9
  circos -conf $ProgDir/Hg_vs_R9_circos.conf -outputdir $OutDir
  mv $OutDir/circos.png $OutDir/Hg_vs_R9_circos.png
  mv $OutDir/circos.svg $OutDir/Hg_vs_R9_circos.svg
  ls $OutDir/Hg_vs_R9_circos.png
```