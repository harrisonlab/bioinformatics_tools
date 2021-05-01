# EPI2ME WIMP 

```bash
srun --partition himem --mem 10G --cpus-per-task 10 --pty bash

cat  /archives/2020_niabemr_nanopore/FOC-test2_S4/20200928_1648_X1_FAO37117_0347a639/fastq_pass/barcode09/*.fastq | gzip > barcode09_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test2_S4/20200928_1648_X1_FAO37117_0347a639/fastq_pass/barcode10/*.fastq | gzip > barcode10_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test2_S4/20200928_1648_X1_FAO37117_0347a639/fastq_pass/barcode08/*.fastq | gzip > barcode08_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test2_S4/20200928_1648_X1_FAO37117_0347a639/fastq_pass/barcode07/*.fastq | gzip > barcode07_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test2_S4/20200928_1648_X1_FAO37117_0347a639/fastq_pass/barcode06/*.fastq | gzip > barcode06_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test2_S4/20200928_1648_X1_FAO37117_0347a639/fastq_pass/barcode05/*.fastq | gzip > barcode05_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test2_S4/20200928_1648_X1_FAO37117_0347a639/fastq_pass/barcode04/*.fastq | gzip > barcode04_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test2_S4/20200928_1648_X1_FAO37117_0347a639/fastq_pass/barcode03/*.fastq | gzip > barcode03_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test2_S4/20200928_1648_X1_FAO37117_0347a639/fastq_pass/barcode02/*.fastq | gzip > barcode02_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test2_S4/20200928_1648_X1_FAO37117_0347a639/fastq_pass/barcode01/*.fastq | gzip > barcode01_pass.fastq.gz

cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/barcode01/*.fastq | gzip > barcode01_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/barcode02/*.fastq | gzip > barcode02_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/barcode03/*.fastq | gzip > barcode03_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/barcode04/*.fastq | gzip > barcode04_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/barcode05/*.fastq | gzip > barcode05_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/barcode08/*.fastq | gzip > barcode08_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/barcode09/*.fastq | gzip > barcode09_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/barcode10/*.fastq | gzip > barcode10_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/barcode11/*.fastq | gzip > barcode11_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/barcode12/*.fastq | gzip > barcode12_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/FOC-test1/20200914_1602_X1_FAN52950_a37c17d4/fastq_pass/unclassified/*.fastq | gzip > unclassified.fastq.gz
```

## minimap2

Reads were aligned to Fus2 and FoL genomes

```bash
    for Reference in $(ls ../../../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa); do
        Read=barcode07_pass.fastq.gz
        OutDir=genome_alignment/minimap2/barcode07/vs_Fus2
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
        sbatch $ProgDir/minimap2.sh $Reference $Read $OutDir
    done

    for Reference in $(ls FoL/GCF_000149955.1_ASM14995v2_genomic.fna); do
        Read=unclassified_pass.fastq.gz
        OutDir=genome_alignment/minimap2/unclassified/vs_FoL
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
        sbatch $ProgDir/minimap2.sh $Reference $Read $OutDir
    done
```

```
# Number of reads
echo $(zcat barcode09_pass.fastq.gz|wc -l)/4|bc
391419
echo $(zcat barcode10_pass.fastq.gz|wc -l)/4|bc
609726

# -F 260 options for primary aligned mapped reads
samtools view -c -F 260 minimap2/Barcode09/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam
280
samtools view -c -F 260 vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam
245

samtools view -c -F 260 minimap2/Barcode10/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam
1833
samtools view -c -F 260 vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam
1798

samtools view -c -F 260 minimap2/unclassified/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam
81
samtools view -c -F 260 vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam
70

samtools view -c -F 260 minimap2/barcode08/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam
162
samtools view -c -F 260 vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam
132

samtools view -c -F 260 minimap2/barcode07/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam
147
samtools view -c -F 260 vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam
103

samtools view -c -F 260 minimap2/barcode06/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam
417
samtools view -c -F 260 vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam
202

samtools view -c -F 260 minimap2/barcode05/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam
169
samtools view -c -F 260 vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam
140
```





samtools view -b -F 4 vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam > mapped.bam
bedtools bamtofastq -i mapped.bam > mapped.fq

samtools view -c -F 260


for Reference in $(ls ncbi-genomes-2020-10-01/GCF_000143535.2_ASM14353v4_genomic.fna); do
Read=Barcode10_pass.fastq.gz
OutDir=genome_alignment/minimap2/barcode10/vs_botrytis
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/minimap2.sh $Reference $Read $OutDir
done

for Reference in $(ls ncbi-genomes-2020-10-01/GCF_000143535.2_ASM14353v4_genomic.fna); do
Read=barcode10_pass.fastq.gz
OutDir=genome_alignment/minimap2/barcode10/vs_botrytis
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/minimap2.sh $Reference $Read $OutDir
done

for Reference in $(ls FoL/GCF_000149955.1_ASM14995v2_genomic.fna); do
Read=barcode10_pass.fastq.gz
OutDir=genome_alignment/minimap2/barcode10/vs_FoL
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/minimap2.sh $Reference $Read $OutDir
done

samtools view -c -F 260 minimap2/barcode10/vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam
1798
samtools view -c -F 260 minimap2/barcode10/vs_botrytis/GCF_000143535.2_ASM14353v4_genomic.fna_aligned_sorted.bam 
79


java -jar /scratch/software/picard-tools-1.89/picard-1.89.jar CompareSAMs \
minimap2/Barcode10/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam \
minimap2/barcode10/vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam \
O=comparison.tsv

samtools view -b -F 260 minimap2/Barcode10/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam | cut -f1 | LC_ALL=C sort | uniq > mapped_reads1_sorted.bam
samtools view -b -F 260 minimap2/barcode10/vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam | cut -f1 | LC_ALL=C sort | uniq > mapped_reads2_sorted.bam
LC_ALL=C comm -12 mapped_reads1_sorted.bam mapped_reads2_sorted.bam > common_seqs.bam



samtools view -F 260 minimap2/Barcode10/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam | cut -f1 | LC_ALL=C sort | uniq > mapped_reads1_sorted.txt
samtools view -F 260 minimap2/barcode10/vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam | cut -f1 | LC_ALL=C sort | uniq > mapped_reads2_sorted.txt
LC_ALL=C comm -12 mapped_reads1_sorted.txt mapped_reads2_sorted.txt > common_seqs.txt
LC_ALL=C comm -23 mapped_reads1_sorted.txt mapped_reads2_sorted.txt > Fus2_uniq.txt
LC_ALL=C comm -13 mapped_reads1_sorted.txt mapped_reads2_sorted.txt > FoL_uniq.txt


1833
1798



1589 




samtools view -b -F 260 minimap2/Barcode10/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam |  cut -f1 | LC_ALL=C sort | uniq  > mapped.bam


samtools view -c -F 260 vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam | cut -f1 | LC_ALL=C sort | uniq > vs_Fus2/mapped_reads1_sorted.txt
samtools view -F 260 vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam | cut -f1 | LC_ALL=C sort | uniq > vs_FoL/mapped_reads2_sorted.txt
LC_ALL=C comm -12 vs_Fus2/mapped_reads1_sorted.txt vs_FoL/mapped_reads2_sorted.txt > common_seqs.txt
LC_ALL=C comm -23 vs_Fus2/mapped_reads1_sorted.txt vs_FoL/mapped_reads2_sorted.txt > Fus2_uniq.txt
LC_ALL=C comm -13 vs_Fus2/mapped_reads1_sorted.txt vs_FoL/mapped_reads2_sorted.txt > FoL_uniq.txt


centrifuge -p 4 -x /data/scratch/gomeza/prog/centrifuge/fungi -U /projects/neonectria_ditissima/gomez_WD/WIMP/barcode10_pass.fastq.gz --min-hitlen 100 --phred33 --report-file centrifuge_report.tsv -S centrifuge_results.txt 

awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' /projects/neonectria_ditissima/gomez_WD/WIMP/barcode10_pass.fastq.gz

centrifuge-kreport -x  /data/scratch/gomeza/prog/centrifuge_v3/fungi centrifuge_results.txt > centrifuge_krakened.txt

echo $(zcat yourfile.fastq.gz|wc -l)/4|bc

```bash
for Reads in $(ls test1/*.fastq.gz); do
echo $Reads
Run=$(echo $Reads | rev | cut -f2 -d '/' | rev)
Sample=$(echo $Reads | rev | cut -f1 -d '/' | rev)
Database=fungi
OutDir=analysis/Centrifuge/$Run/$Sample
echo $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
sbatch $ProgDir/centrifuge.sh $Database $OutDir $Reads
done
```

```bash
for Read1 in $(ls Charlotte/star_aligmentUnmapped.out.mate1); do
Read2=$(echo $Read1 | sed 's/mate1/mate2/g')
echo $Read1
echo $Read2
Sample=$(echo $Read1 | rev | cut -f2 -d '/' | rev)
Database=plantvirus
OutDir=analysis/Centrifuge_vtest/$Sample
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
sbatch $ProgDir/centrifuge_v2.sh $Database $OutDir $Read1 $Read2
done
```
```

srun --partition short --mem-per-cpu 6G --cpus-per-task 60 --pty bash



```bash
srun --partition himem --mem 10G --cpus-per-task 10 --pty bash

cat  /archives/2020_niabemr_nanopore/S3D6/20201215_1453_X1_FAO27987_9702022c/fastq_pass/*.fastq | gzip > raw_data/S3D6_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-2/20201215_1645_X1_FAO27987_b2861159/fastq_pass/*.fastq | gzip > raw_data/S3D6-2_pass.fastq.gz

cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode01/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode01_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode02/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode02_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode03/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode03_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode04/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode04_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode05/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode05_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode06/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode06_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode07/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode07_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode08/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode08_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode09/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode09_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode10/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode10_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode11/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode11_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/barcode12/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/barcode12_pass.fastq.gz
cat  /archives/2020_niabemr_nanopore/S3D6-BC11_BC12_No_Selection/20201217_1441_X2_FAO27999_57821556/fastq_pass/unclassified/*.fastq | gzip > raw_data/S3D6-BC11_BC12_No_Selection/unclassified_pass.fastq.gz


cat  /archives/2020_niabemr_nanopore/S3_D6_readfish1/20201216_1220_X1_FAO27987_4cc4da2e/fastq_pass/*.fastq | gzip > raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
```
```bash
# Run fastqc
RawData=$(ls raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz)
OutDir=qc_data/S3_D6_readfish1
echo $RawData;
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
sbatch $ProgDir/fastqc2.sh $RawData $OutDir

```

## minimap2

Reads were aligned to Fus2 and FoL genomes

```bash
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa); do
Read=raw_data/S3D6-BC11_BC12_No_Selection/barcode12_pass.fastq.gz
OutDir=genome_alignment/minimap2/S3D6_NoSelection/barcode12/vs_Fus2
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch -p long $ProgDir/minimap2.sh $Reference $Read $OutDir
done
```


echo $(zcat S3D6_pass.fastq.gz|wc -l)/4|bc
125373

echo $(zcat S3D6-2_pass.fastq.gz|wc -l)/4|bc
11003176

echo $(zcat S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz|wc -l)/4|bc
9764274

echo $(zcat S3D6-BC11_BC12_No_Selection/barcode01_pass.fastq.gz|wc -l)/4|bc
9
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode02_pass.fastq.gz|wc -l)/4|bc
10
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode03_pass.fastq.gz|wc -l)/4|bc
43
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode04_pass.fastq.gz|wc -l)/4|bc
9
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode05_pass.fastq.gz|wc -l)/4|bc
8
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode06_pass.fastq.gz|wc -l)/4|bc
26
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode07_pass.fastq.gz|wc -l)/4|bc
14
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode08_pass.fastq.gz|wc -l)/4|bc
11
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode09_pass.fastq.gz|wc -l)/4|bc
15
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode10_pass.fastq.gz|wc -l)/4|bc
2
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode11_pass.fastq.gz|wc -l)/4|bc
910153
echo $(zcat S3D6-BC11_BC12_No_Selection/barcode12_pass.fastq.gz|wc -l)/4|bc
1322587
echo $(zcat S3D6-BC11_BC12_No_Selection/unclassified_pass.fastq.gz|wc -l)/4|bc
29044



# -F 260 options for primary aligned mapped reads
samtools view -c genome_alignment/minimap2/S3D6/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
125637
#mapped reads
samtools view -c -F 4 genome_alignment/minimap2/S3D6/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
632
samtools view -c -f 4 genome_alignment/minimap2/S3D6/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
125005
samtools view -c -F 260 genome_alignment/minimap2/S3D6/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
421
samtools view -F 0x904 -c genome_alignment/minimap2/S3D6/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
368

samtools view -c -F 4 genome_alignment/minimap2/S3D6-2/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
33320
samtools view -c -F 260 genome_alignment/minimap2/S3D6-2/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
23911

samtools view -c -F 4 genome_alignment/minimap2/S3_D6_readfish/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
29483
samtools view -c -F 260 genome_alignment/minimap2/S3_D6_readfish/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
20298

samtools view -c -F 4 genome_alignment/minimap2/S3D6_NoSelection/barcode11/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
4528
samtools view -c -F 260 genome_alignment/minimap2/S3D6_NoSelection/barcode11/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
2851

samtools view -c -F 4 genome_alignment/minimap2/S3D6_NoSelection/barcode12/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam

samtools view -c -F 260 genome_alignment/minimap2/S3D6_NoSelection/barcode12/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
4228

samtools view -c -F 260 genome_alignment/minimap2/S3D6_NoSelection/unclassified/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam
103


samtools view -c -F 4 genome_alignment/minimap2/S3D6/vs_Fus2/sorted.bam 


bedtools genomecov -ibam genome_alignment/minimap2/S3D6/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned.bam -bga > genome_alignment/minimap2/S3D6/vs_Fus2/output_aln.bam
awk '{s+=$4}END{print s}' output_aln.bam



bedtools genomecov -ibam genome_alignment/minimap2/S3D6/vs_Fus2/sorted.bam -bga > genome_alignment/minimap2/S3D6/vs_Fus2/output_aln.bam
awk '{s+=$4}END{print s}' genome_alignment/minimap2/S3D6/vs_Fus2/output_aln.bam
39287



bedtools genomecov -ibam genome_alignment/minimap2/S3D6-2/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam -bga > genome_alignment/minimap2/S3D6-2/vs_Fus2/output_aln.bam
awk '{s+=$4}END{print s}' genome_alignment/minimap2/S3D6-2/vs_Fus2/output_aln.bam


bedtools genomecov -ibam genome_alignment/minimap2/S3_D6_readfish/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam -bga > genome_alignment/minimap2/S3_D6_readfish/vs_Fus2/output_aln.bam
awk '{s+=$4}END{print s}' genome_alignment/minimap2/S3_D6_readfish/vs_Fus2/output_aln.bam

bedtools genomecov -ibam genome_alignment/minimap2/S3D6_NoSelection/barcode11/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam -bga > genome_alignment/minimap2/S3D6_NoSelection/barcode11/vs_Fus2/output_aln.bam
awk '{s+=$4}END{print s}' genome_alignment/minimap2/S3D6_NoSelection/barcode11/vs_Fus2/output_aln.bam

bedtools genomecov -ibam genome_alignment/minimap2/S3D6_NoSelection/barcode12/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam -bga > genome_alignment/minimap2/S3D6_NoSelection/barcode12/vs_Fus2/output_aln.bam
awk '{s+=$4}END{print s}' genome_alignment/minimap2/S3D6_NoSelection/barcode12/vs_Fus2/output_aln.bam

bedtools genomecov -ibam genome_alignment/minimap2/S3D6_NoSelection/barcode12/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam -bg > output_aln.bg
awk '{ s+=($3-$2) } END { print s }' output_aln.bg

# Estimate coverage long read data
for RawData in $(ls -d raw_data/S3D6-2_pass.fastq.gz); do
echo $RawData
GenomeSize=53 #Estimated genome size
OutDir=$(dirname $RawData)
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
sbatch $ProgDir/count_nucl_single.sh $RawData $GenomeSize $OutDir
done
```

cat raw_data/S3D6-BC11_BC12_No_Selection/barcode12_pass.fastq.gz | gunzip -fc > 05.fastq
/data/scratch/gomeza/prog/count_nucl.pl -i 05.fastq -g 53 > 05.cov



## minimap2

Reads were aligned to Fus2 and FoL genomes

```bash
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_mathioli/Stocks4/filtered_contigs/Stocks4_contigs_unmasked.fa); do
Read=raw_data/S3D6_pass.fastq.gz
OutDir=genome_alignment/minimap2/S3D6/vs_FoM
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch -p long $ProgDir/minimap2.sh $Reference $Read $OutDir
done
```

```bash
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl_repmask/4287_chromosomal_contigs_unmasked.fa); do
Read=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
OutDir=genome_alignment/minimap2/S3_D6_readfish/vs_FoL
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch -p long $ProgDir/minimap2.sh $Reference $Read $OutDir
done
```

samtools view -c -F 4 genome_alignment/minimap2/S3D6/vs_FoL/4287_chromosomal_contigs_unmasked.fa_aligned_sorted.bam

samtools view -c -F 260 genome_alignment/minimap2/S3D6/vs_FoL/4287_chromosomal_contigs_unmasked.fa_aligned_sorted.bam

samtools view -c -F 4 genome_alignment/minimap2/S3D6-2/vs_FoL/4287_chromosomal_contigs_unmasked.fa_aligned_sorted.bam

samtools view -c -F 260 genome_alignment/minimap2/S3D6-2/vs_FoL/4287_chromosomal_contigs_unmasked.fa_aligned_sorted.bam

samtools view -c -F 4 genome_alignment/minimap2/S3_D6_readfish/vs_FoL/4287_chromosomal_contigs_unmasked.fa_aligned_sorted.bam

samtools view -c -F 260 genome_alignment/minimap2/S3_D6_readfish/vs_FoL/4287_chromosomal_contigs_unmasked.fa_aligned_sorted.bam


python /data/scratch/gomeza/prog/LongQC/longQC.py sampleqc -x ont-rapid -o out_dir S3D6-BC11_BC12_No_Selection/barcode01_pass.fastq.gz

NanoFilt --logfile log.txt -l 90 -q 10 04.fastq >  highQuality-reads.fastq