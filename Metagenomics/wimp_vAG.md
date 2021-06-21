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


## Fusarium genomes

```bash
# F.oxysporum
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum/fo47/ncbi_edits_repmask/fo47_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum/fo47/broad_repmask/fo47_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum/fo47_tgac_filtered/ncbi_edits_repmask/fo47_tgac_filtered_contigs_unmasked.fa

# F.oxysporum fsp cepae
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_f.sp.cepae/55/ncbi_edits_repmask/55_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_f.sp.cepae/55_tgac_filtered/ncbi_edits_repmask/55_tgac_filtered_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/125_ncbi/ncbi_submission/125_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/A13_ncbi/ncbi_submission/A13_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/A23_ncbi/ncbi_submission/A23_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/A28_ncbi/ncbi_submission/A28_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/CB3_ncbi/ncbi_submission/CB3_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/PG_ncbi/ncbi_submission/PG_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa

../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae_old/A1-2/filtered_contigs_repmask/A1-2_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae_old/D2/filtered_contigs_repmask/D2_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae_old/HB6/filtered_contigs_repmask/HB6_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae_old/HB17/filtered_contigs_repmask/HB17_contigs_unmasked.fa

# F.oxysporum fsp fragariae
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_fragariae/15-074/filtered_contigs/15-074_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_fragariae/Straw465/filtered_contigs/Straw465_contigs_unmasked.fa

# F.oxysporum fsp lactucae
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lactucae/AJ516/filtered_contigs/AJ516_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lactucae/AJ516/filtered_ncbi/AJ516_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lactucae/AJ520/filtered_contigs/AJ520_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lactucae/AJ520/filtered_ncbi/AJ520_contigs_unmasked.fa

# F.oxysporum fsp lycopersici
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl_repmask/4287_chromosomal_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa

# F.oxysporum fsp mathioli
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_mathioli/Stocks4/filtered_contigs/Stocks4_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_mathioli/Stocks4/filtered_ncbi/Stocks4_contigs_unmasked.fa

# F.oxysporum narcissi
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_narcissi/FON129/filtered_contigs/FON129_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_narcissi/FON_63/filtered_contigs/FON_63_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_narcissi/FON77/filtered_contigs/FON77_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_narcissi/FON81/filtered_contigs/FON81_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_narcissi/FON89/filtered_contigs/FON89_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_narcissi/N139_ncbi/ncbi_submission/N139_contigs_unmasked.fa

# F.oxysporum fsp pisi
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_pisi/F81/filtered_contigs/F81_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_pisi/FOP1/filtered_contigs_repmask/FOP1_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_pisi/FOP1-EMR/filtered_contigs/FOP1-EMR_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_pisi/FOP2/filtered_contigs_repmask/FOP2_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_pisi/FOP5/filtered_contigs_repmask/FOP5_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_pisi/PG18/filtered_contigs_repmask/PG18_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_pisi/PG3/filtered_contigs_repmask/PG3_contigs_unmasked.fa
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_pisi/PG3/filtered_contigs_repmask/PG3_contigs_unmasked.fa

# F.oxysporum fsp statice
../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_statice/Stat10/filtered_contigs/Stat10_contigs_unmasked.fa
```


## minimap2


```bash

# Fo
for Strain in fo47 fo47_tgac_filtered; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/ncb*/*_contigs_unmasked.fa); do
echo $Reference
Read=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
R1=$(echo $Read | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R1
Read2=raw_data/S3D6_pass.fastq.gz
R2=$(echo $Read2 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R2
Read3=raw_data/S3D6-2_pass.fastq.gz
R3=$(echo $Read3 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/Fo/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
OutDir2=genome_alignment/minimap2/$R2/Fo/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read2 $OutDir2
OutDir3=genome_alignment/minimap2/$R3/Fo/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read3 $OutDir3
for bar in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
Read4=raw_data/S3D6-BC11_BC12_No_Selection/"$bar"_pass.fastq.gz
R4=$(echo $Read4 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R4
OutDir4=genome_alignment/minimap2/S3D6_NoSelection/$R4/Fo/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read4 $OutDir4
done
done
done

# Fo cepae
for Strain in 55 55_tgac_filtered 125_ncbi A13_ncbi A23_ncbi A28_ncbi CB3_ncbi PG_ncbi Fus2_canu_new; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
R1=$(echo $Read | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R1
Read2=raw_data/S3D6_pass.fastq.gz
R2=$(echo $Read2 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R2
Read3=raw_data/S3D6-2_pass.fastq.gz
R3=$(echo $Read3 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoC/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
OutDir2=genome_alignment/minimap2/$R2/FoF/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read2 $OutDir2
OutDir3=genome_alignment/minimap2/$R3/FoF/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read3 $OutDir3
for bar in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
Read4=raw_data/S3D6-BC11_BC12_No_Selection/"$bar"_pass.fastq.gz
R4=$(echo $Read4 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R4
OutDir4=genome_alignment/minimap2/$R4/FoF/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read4 $OutDir4
done
done
done



# F.oxysporum fsp fragariae
for Strain in 15-074 Straw465; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/filtered_contigs/*_contigs_unmasked.fa); do
echo $Reference
Read=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
R1=$(echo $Read | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R1
Read2=raw_data/S3D6_pass.fastq.gz
R2=$(echo $Read2 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R2
Read3=raw_data/S3D6-2_pass.fastq.gz
R3=$(echo $Read3 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoF/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
OutDir2=genome_alignment/minimap2/$R2/FoF/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read2 $OutDir2
OutDir3=genome_alignment/minimap2/$R3/FoF/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read3 $OutDir3
for bar in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
Read4=raw_data/S3D6-BC11_BC12_No_Selection/"$bar"_pass.fastq.gz
R4=$(echo $Read4 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R4
OutDir4=genome_alignment/minimap2/$R4/FoF/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read4 $OutDir4
done
done
done

# Fo lactucae
for Strain in AJ516 AJ520; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/filtered_contigs/*_contigs_unmasked.fa); do
echo $Reference
Read=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
R1=$(echo $Read | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R1
Read2=raw_data/S3D6_pass.fastq.gz
R2=$(echo $Read2 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R2
Read3=raw_data/S3D6-2_pass.fastq.gz
R3=$(echo $Read3 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoL/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
OutDir2=genome_alignment/minimap2/$R2/FoL/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read2 $OutDir2
OutDir3=genome_alignment/minimap2/$R3/FoL/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read3 $OutDir3
for bar in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
Read4=raw_data/S3D6-BC11_BC12_No_Selection/"$bar"_pass.fastq.gz
R4=$(echo $Read4 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R4
OutDir4=genome_alignment/minimap2/$R4/FoL/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read4 $OutDir4
done
done
done


# Fo lycopersici
for Strain in 4287_chromosomal 4287_v2; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
R1=$(echo $Read | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R1
Read2=raw_data/S3D6_pass.fastq.gz
R2=$(echo $Read2 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R2
Read3=raw_data/S3D6-2_pass.fastq.gz
R3=$(echo $Read3 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoLy/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
OutDir2=genome_alignment/minimap2/$R2/FoLy/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read2 $OutDir2
OutDir3=genome_alignment/minimap2/$R3/FoLy/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read3 $OutDir3
for bar in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
Read4=raw_data/S3D6-BC11_BC12_No_Selection/"$bar"_pass.fastq.gz
R4=$(echo $Read4 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R4
OutDir4=genome_alignment/minimap2/S3D6_NoSelection/$R4/FoLy/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read4 $OutDir4
done
done
done


# Fo mathioli
for Strain in Stocks4; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/filtered_contigs/*_contigs_unmasked.fa); do
echo $Reference
Read=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
R1=$(echo $Read | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R1
Read2=raw_data/S3D6_pass.fastq.gz
R2=$(echo $Read2 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R2
Read3=raw_data/S3D6-2_pass.fastq.gz
R3=$(echo $Read3 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoM/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
OutDir2=genome_alignment/minimap2/$R2/FoM/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read2 $OutDir2
OutDir3=genome_alignment/minimap2/$R3/FoM/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read3 $OutDir3
for bar in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
Read4=raw_data/S3D6-BC11_BC12_No_Selection/"$bar"_pass.fastq.gz
R4=$(echo $Read4 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R4
OutDir4=genome_alignment/minimap2/S3D6_NoSelection/$R4/FoM/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read4 $OutDir4
done
done
done

# F.oxysporum narcissi
for Strain in FON129 FON_63 FON77 FON81 FON89 N139_ncbi; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
R1=$(echo $Read | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R1
Read2=raw_data/S3D6_pass.fastq.gz
R2=$(echo $Read2 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R2
Read3=raw_data/S3D6-2_pass.fastq.gz
R3=$(echo $Read3 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoN/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
OutDir2=genome_alignment/minimap2/$R2/FoN/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read2 $OutDir2
OutDir3=genome_alignment/minimap2/$R3/FoN/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read3 $OutDir3
for bar in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
Read4=raw_data/S3D6-BC11_BC12_No_Selection/"$bar"_pass.fastq.gz
R4=$(echo $Read4 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R4
OutDir4=genome_alignment/minimap2/S3D6_NoSelection/$R4/FoN/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read4 $OutDir4
done
done
done


# F.oxysporum fsp pisi
for Strain in F81 FOP1 FOP1-EMR FOP2 FOP5 PG18 PG3; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
R1=$(echo $Read | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R1
Read2=raw_data/S3D6_pass.fastq.gz
R2=$(echo $Read2 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R2
Read3=raw_data/S3D6-2_pass.fastq.gz
R3=$(echo $Read3 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoP/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
OutDir2=genome_alignment/minimap2/$R2/FoP/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read2 $OutDir2
OutDir3=genome_alignment/minimap2/$R3/FoP/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read3 $OutDir3
for bar in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
Read4=raw_data/S3D6-BC11_BC12_No_Selection/"$bar"_pass.fastq.gz
R4=$(echo $Read4 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R4
OutDir4=genome_alignment/minimap2/S3D6_NoSelection/$R4/FoP/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read4 $OutDir4
done
done
done


# F.oxysporum fsp statice
for Strain in Stat10; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz
R1=$(echo $Read | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R1
Read2=raw_data/S3D6_pass.fastq.gz
R2=$(echo $Read2 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R2
Read3=raw_data/S3D6-2_pass.fastq.gz
R3=$(echo $Read3 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R3
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoS/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
OutDir2=genome_alignment/minimap2/$R2/FoS/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read2 $OutDir2
OutDir3=genome_alignment/minimap2/$R3/FoS/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read3 $OutDir3
for bar in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
Read4=raw_data/S3D6-BC11_BC12_No_Selection/"$bar"_pass.fastq.gz
R4=$(echo $Read4 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R4
OutDir4=genome_alignment/minimap2/S3D6_NoSelection/$R4/FoS/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read4 $OutDir4
done
done
done
```

for Reads in $(ls minimap2/S3_D6_readfish/*/*/*_contigs_unmasked.fa_aligned_sorted.bam ); do
Org=$(echo $Reads | rev | cut -f3 -d '/' | rev)
echo $Org
vs=$(echo $Reads | rev | cut -f2 -d '/' | rev)
echo $vs
samtools view -F 260 $Reads | cut -f1 | LC_ALL=C sort | uniq > minimap2/S3_D6_readfish/$Org/$vs/mapped_"$vs"_sorted.txt
done


for Reads in $(ls minimap2/S3D6-2/*/*/*_contigs_unmasked.fa_aligned_sorted.bam ); do
Org=$(echo $Reads | rev | cut -f3 -d '/' | rev)
echo $Org
vs=$(echo $Reads | rev | cut -f2 -d '/' | rev)
echo $vs
samtools view -F 260 $Reads | cut -f1 | LC_ALL=C sort | uniq > minimap2/S3D6-2/$Org/$vs/mapped_"$vs"_sorted.txt
done


for Reads in $(ls minimap2/S3D6_NoSelection/barcode11/*/*/*_contigs_unmasked.fa_aligned_sorted.bam ); do
Org=$(echo $Reads | rev | cut -f3 -d '/' | rev)
echo $Org
vs=$(echo $Reads | rev | cut -f2 -d '/' | rev)
echo $vs
samtools view -F 260 $Reads | cut -f1 | LC_ALL=C sort | uniq > minimap2/S3D6_NoSelection/barcode11/$Org/$vs/mapped_"$vs"_sorted.txt
done

for Reads in $(ls minimap2/S3D6_NoSelection/barcode12/*/*/*_contigs_unmasked.fa_aligned_sorted.bam ); do
Org=$(echo $Reads | rev | cut -f3 -d '/' | rev)
echo $Org
vs=$(echo $Reads | rev | cut -f2 -d '/' | rev)
echo $vs
samtools view -F 260 $Reads | cut -f1 | LC_ALL=C sort | uniq > minimap2/S3D6_NoSelection/barcode12/$Org/$vs/mapped_"$vs"_sorted.txt
done




samtools view -F 260 minimap2/Barcode10/vs_Fus2/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam | cut -f1 | LC_ALL=C sort | uniq > mapped_reads1_sorted.txt
samtools view -F 260 minimap2/barcode10/vs_FoL/GCF_000149955.1_ASM14995v2_genomic.fna_aligned_sorted.bam | cut -f1 | LC_ALL=C sort | uniq > mapped_reads2_sorted.txt
LC_ALL=C comm -12 mapped_reads1_sorted.txt mapped_reads2_sorted.txt > common_seqs.txt
LC_ALL=C comm -23 mapped_reads1_sorted.txt mapped_reads2_sorted.txt > Fus2_uniq.txt
LC_ALL=C comm -13 mapped_reads1_sorted.txt mapped_reads2_sorted.txt > FoL_uniq.txt


LC_ALL=C comm -12 mapped_reads1_sorted.txt mapped_reads2_sorted.txt > common_seqs.txt


LC_ALL=C comm -23 all_S3D6-2_FoC.txt all_S3D6-2_FoF.txt  > FOC_uniq1.txt
LC_ALL=C comm -23 all_S3D6-2_FoC.txt all_S3D6-2_FoL.txt  > FOC_uniq2.txt
LC_ALL=C comm -23 all_S3D6-2_FoC.txt all_S3D6-2_FoLy.txt  > FOC_uniq3.txt
LC_ALL=C comm -23 all_S3D6-2_FoC.txt all_S3D6-2_FoM.txt  > FOC_uniq4.txt
LC_ALL=C comm -23 all_S3D6-2_FoC.txt all_S3D6-2_FoN.txt  > FOC_uniq5.txt
LC_ALL=C comm -23 all_S3D6-2_FoC.txt all_S3D6-2_FoP.txt  > FOC_uniq6.txt
LC_ALL=C comm -23 all_S3D6-2_FoC.txt all_S3D6-2_FoS.txt  > FOC_uniq7.txt
LC_ALL=C comm -23 all_S3D6-2_FoC.txt all_S3D6-2_Fo.txt > FOC_uniq8.txt

cat FOC* | sort | uniq > FOC_unique_all.txt


LC_ALL=C comm -23 all_barcode12_FoC.txt all_barcode12_FoF.txt  > FOC_uniq1.txt
LC_ALL=C comm -23 all_barcode12_FoC.txt all_barcode12_FoL.txt  > FOC_uniq2.txt
LC_ALL=C comm -23 all_barcode12_FoC.txt all_barcode12_FoLy.txt  > FOC_uniq3.txt
LC_ALL=C comm -23 all_barcode12_FoC.txt all_barcode12_FoM.txt  > FOC_uniq4.txt
LC_ALL=C comm -23 all_barcode12_FoC.txt all_barcode12_FoN.txt  > FOC_uniq5.txt
LC_ALL=C comm -23 all_barcode12_FoC.txt all_barcode12_FoP.txt  > FOC_uniq6.txt
LC_ALL=C comm -23 all_barcode12_FoC.txt all_barcode12_FoS.txt  > FOC_uniq7.txt
LC_ALL=C comm -23 all_barcode12_FoC.txt all_barcode12_Fo.txt > FOC_uniq8.txt


LC_ALL=C comm -13 mapped_reads1_sorted.txt mapped_reads2_sorted.txt > FoL_uniq.txt




cat FoC/*/mapped_*.txt | uniq > all_S3_D6_readfish_FoC.txt 
cat S3_D6_readfish/*/*/mapped_*.txt | uniq > S3_D6_readfish/all_S3_D6_readfish.txt 
cat S3_D6_readfish/*/*/mapped_*.txt | uniq > S3_D6_readfish/all_S3_D6_readfish.txt 
cat S3_D6_readfish/*/*/mapped_*.txt | uniq > S3_D6_readfish/all_S3_D6_readfish.txt 


```bash
for D1 in $(ls *contigs_unmasked.fa) do
R1=$(echo $D1 | rev | cut -f1 -d '/' | rev | sed 's/_contigs_unmasked.fa//g')
cat $D1 | sed 's/contig_/Fo_"$R1"_contig_/g' > "$R1"_contigs_unmasked_renamed.fa
done

less Fus2_canu_contigs_unmasked.fa | sed 's/contig_/Fo_cepae_Fus2_contig_/g' | sed 's/_pilon//g' > Fus2_contigs_unmasked_renamed.fa
less Straw465_renamed.fasta  | sed 's/contig_/Fo_fragariae_Straw465_contig_/g'  > Straw465_contigs_unmasked_renamed.fa
cat 


# read histogram

RawData=$(ls raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz)
echo $RawData;
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
sbatch -p long $ProgDir/fastqc.sh $RawData


readlength.sh in=raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz out=histogram.txt

## extract reads!!!!

/scratch/software/seqtk/seqtk subseq raw_data/S3_D6_readfish1/S3D6_readfish.fastq /archives/2020_niabemr_nanopore/S3_D6_readfish1/20201216_1220_X1_FAO27987_4cc4da2e/unblocked_read_ids.txt > out.fq


echo $(cat out.fq |wc -l)/4|bc

awk '{if(NR%4==2) print length($1)}' out.fq | sort -n | uniq -c > read_length2.txt

reads<-read.csv(file="read_length2.txt", sep="", header=FALSE)
plot (reads$V2,reads$V1,type="l",xlab="read length",ylab="occurences",col="blue")

/scratch/software/seqkit/seqkit seq -m 1000 raw_data/S3_D6_readfish1/S3D6_readfish.fastq > filtered.fq
echo $(cat filtered.fq|wc -l)/4|bc


for Read1 in $(ls filtered.fq.gz ); do
Out=$(echo $Read1 | rev | cut -f1 -d '/' | rev)
echo $Out
Database=Foxysporum
OutDir=analysis/Kraken/$Out
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
sbatch $ProgDir/kraken_long_reads.sh $Read1 $Database $OutDir
done


awk '{if(NR%4==2) print length($1)}' classifiedseqs.fq | sort -n | uniq -c > classified_read_length.txt
awk '{if(NR%4==2) print length($1)}' unclassifiedseqs.fq | sort -n | uniq -c > unclassified_read_length.txt


reads<-read.csv(file="classified_read_length.txt", sep="", header=FALSE)
plot (reads$V2,reads$V1,type="l",xlab="read length",ylab="occurences",col="blue")

reads<-read.csv(file="unclassified_read_length.txt", sep="", header=FALSE)
plot (reads$V2,reads$V1,type="l",xlab="read length",ylab="occurences",col="blue")


/scratch/software/seqtk/seqtk subseq classifiedseqs.fq /archives/2020_niabemr_nanopore/S3_D6_readfish1/20201216_1220_X1_FAO27987_4cc4da2e/unblocked_read_ids.txt > classified_unblocked.fq
/scratch/software/seqtk/seqtk subseq unclassifiedseqs.fq /archives/2020_niabemr_nanopore/S3_D6_readfish1/20201216_1220_X1_FAO27987_4cc4da2e/unblocked_read_ids.txt > unclassified_unblocked.fq


echo $(cat classified_unblocked.fq|wc -l)/4|bc
8
echo $(cat unclassified_unblocked.fq|wc -l)/4|bc
628

echo $(cat classifiedseqs.fq|wc -l)/4|bc
5222
echo $(cat unclassifiedseqs.fq|wc -l)/4|bc
3377



for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=filtered.fq.gz
R1=$(echo $Read | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $R1
OutDir=genome_alignment/minimap2/$R1/FoS/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
done


```bash

# Fo
for Strain in fo47 fo47_tgac_filtered; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/ncb*/*_contigs_unmasked.fa); do
echo $Reference
Read=filtered.fq.gz
R1=S3D6_readfish_filtered
echo $R1
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/Fo/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
done
done

# Fo cepae
for Strain in 55 55_tgac_filtered 125_ncbi A13_ncbi A23_ncbi A28_ncbi CB3_ncbi PG_ncbi Fus2_canu_new; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/ed*/*_contigs_unmasked.fa); do
echo $Reference
Read=filtered.fq.gz
R1=S3D6_readfish_filtered
echo $R1
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoC/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
done
done

# F.oxysporum fsp fragariae
for Strain in 15-074 Straw465; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=filtered.fq.gz
R1=S3D6_readfish_filtered
echo $R1
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoF/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
done
done

# Fo lactucae
for Strain in AJ516 AJ520; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=filtered.fq.gz
R1=S3D6_readfish_filtered
echo $R1
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoL/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
done
done

# Fo lycopersici
for Strain in 4287_chromosomal 4287_v2; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=filtered.fq.gz
R1=S3D6_readfish_filtered
echo $R1
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoLy/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
done
done


# Fo mathioli
for Strain in Stocks4; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=filtered.fq.gz
R1=S3D6_readfish_filtered
echo $R1
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoM/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
done
done

# F.oxysporum narcissi
for Strain in FON129 FON_63 FON77 FON81 FON89 N139_ncbi; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=filtered.fq.gz
R1=S3D6_readfish_filtered
echo $R1
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoN/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
done
done


# F.oxysporum fsp pisi
for Strain in F81 FOP1 FOP1-EMR FOP2 FOP5 PG18 PG3; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=filtered.fq.gz
R1=S3D6_readfish_filtered
echo $R1
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoP/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
done
done


# F.oxysporum fsp statice
for Strain in Stat10; do
for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/*/$Strain/*/*_contigs_unmasked.fa); do
echo $Reference
Read=filtered.fq.gz
R1=S3D6_readfish_filtered
echo $R1
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
OutDir=genome_alignment/minimap2/$R1/FoS/vs_"$Strain"
sbatch -p short $ProgDir/minimap2.sh $Reference $Read $OutDir
done
done

for Reads in $(ls genome_alignment/minimap2/S3D6_readfish_filtered/*/*/*_contigs_unmasked.fa_aligned_sorted.bam ); do
Org=$(echo $Reads | rev | cut -f3 -d '/' | rev)
echo $Org
vs=$(echo $Reads | rev | cut -f2 -d '/' | rev)
echo $vs
samtools view -F 260 $Reads | cut -f1 | LC_ALL=C sort | uniq > genome_alignment/minimap2/S3D6_readfish_filtered/$Org/$vs/mapped_"$vs"_sorted.txt
done

for Reads in $(ls genome_alignment/minimap2/S3D6_readfish_filtered/*/*/mapped_*_sorted.txt); do
Org=$(echo $Reads | rev | cut -f2 -d '/' | rev)
echo $Org
less $Reads | wc -l
done

cd genome_alignment/minimap2/S3D6_readfish_filtered

cat Fo/*/mapped_*_sorted.txt | sort | uniq  > all_Fo_reads.txt
cat FoC/*/mapped_*_sorted.txt | sort | uniq  > all_FoC_reads.txt
cat FoF/*/mapped_*_sorted.txt | sort | uniq  > all_FoF_reads.txt
cat FoL/*/mapped_*_sorted.txt | sort | uniq  > all_FoL_reads.txt
cat FoLy/*/mapped_*_sorted.txt | sort | uniq  > all_FoLy_reads.txt
cat FoM/*/mapped_*_sorted.txt | sort | uniq  > all_FoM_reads.txt
cat FoN/*/mapped_*_sorted.txt | sort | uniq  > all_FoN_reads.txt
cat FoP/*/mapped_*_sorted.txt | sort | uniq  > all_FoP_reads.txt
cat FoS/*/mapped_*_sorted.txt | sort | uniq  > all_FoS_reads.txt

for d in all_Fo*; do
less $d | wc -l
done

5175
5054
5045
4903
4938
5131
5168
4987
4930

LC_ALL=C comm -23 all_FoC_reads.txt all_Fo_reads.txt  > FOC_uniq1.txt
LC_ALL=C comm -23 all_FoC_reads.txt all_FoF_reads.txt  > FOC_uniq2.txt
LC_ALL=C comm -23 all_FoC_reads.txt all_FoL_reads.txt  > FOC_uniq3.txt
LC_ALL=C comm -23 all_FoC_reads.txt all_FoLy_reads.txt  > FOC_uniq4.txt
LC_ALL=C comm -23 all_FoC_reads.txt all_FoM_reads.txt  > FOC_uniq5.txt
LC_ALL=C comm -23 all_FoC_reads.txt all_FoN_reads.txt  > FOC_uniq6.txt
LC_ALL=C comm -23 all_FoC_reads.txt all_FoP_reads.txt  > FOC_uniq7.txt
LC_ALL=C comm -23 all_FoC_reads.txt all_FoS_reads.txt > FOC_uniq8.txt

for d in FOC_uniq* ; do
less $d | wc -l
done

189
121
132
273
239
98
31
246

cat FOC_uniq* 
```

## Rename contigs

```bash
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
# If split or remove contigs is needed, provide FCSreport file by NCBI.
touch tmp.txt
for Assembly in $(ls R2_contigs_unmasked.fa); do
OutDir=$(dirname $Assembly)
$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/R2_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
rm tmp.txt
```

kraken2-build --download-taxonomy --db Foxysporum
kraken2-build --build --db Foxysporum
c
srun --partition long --mem-per-cpu 10G --cpus-per-task 10 --pty bash


```bash

for file in $(ls fusarium_db/FoCepae/*_contigs_unmasked_renamed.fa); do
R1=$(echo $file | rev | cut -f1 -d '/' | rev | sed 's/_renamed.fa//g')
perl /home/gomeza/miniconda3/envs/Metagenomics/bin/util/addTaxonIDToFasta.pl --inputFA $file --outputFA "$R1"_kraken.fa --taxonID 396571
done

for file in $(ls fusarium_db/Fo/*_contigs_unmasked_renamed.fa); do
R1=$(echo $file | rev | cut -f1 -d '/' | rev | sed 's/_renamed.fa//g')
perl /home/gomeza/miniconda3/envs/Metagenomics/bin/util/addTaxonIDToFasta.pl --inputFA $file --outputFA "$R1"_kraken.fa --taxonID 660027
done

for file in $(ls fusarium_db/FoFragariae/*_contigs_unmasked_renamed.fa); do
R1=$(echo $file | rev | cut -f1 -d '/' | rev | sed 's/_renamed.fa//g')
perl /home/gomeza/miniconda3/envs/Metagenomics/bin/util/addTaxonIDToFasta.pl --inputFA $file --outputFA "$R1"_kraken.fa --taxonID 100903
done

for file in $(ls fusarium_db/Folactucae/*_contigs_unmasked_renamed.fa); do
R1=$(echo $file | rev | cut -f1 -d '/' | rev | sed 's/_renamed.fa//g')
perl /home/gomeza/miniconda3/envs/Metagenomics/bin/util/addTaxonIDToFasta.pl --inputFA $file --outputFA "$R1"_kraken.fa --taxonID 299031
done


for file in $(ls fusarium_db/FoLycopersici/*_contigs_unmasked_renamed.fa); do
R1=$(echo $file | rev | cut -f1 -d '/' | rev | sed 's/_renamed.fa//g')
perl /home/gomeza/miniconda3/envs/Metagenomics/bin/util/addTaxonIDToFasta.pl --inputFA $file --outputFA "$R1"_kraken.fa --taxonID 4264281
done

for file in $(ls fusarium_db/FoNarcissi/*_contigs_unmasked_renamed.fa); do
R1=$(echo $file | rev | cut -f1 -d '/' | rev | sed 's/_renamed.fa//g')
perl /home/gomeza/miniconda3/envs/Metagenomics/bin/util/addTaxonIDToFasta.pl --inputFA $file --outputFA "$R1"_kraken.fa --taxonID 451672
done

for file in $(ls fusarium_db/FoPisi/*_contigs_unmasked_renamed.fa); do
R1=$(echo $file | rev | cut -f1 -d '/' | rev | sed 's/_renamed.fa//g')
perl /home/gomeza/miniconda3/envs/Metagenomics/bin/util/addTaxonIDToFasta.pl --inputFA $file --outputFA "$R1"_kraken.fa --taxonID 179143
done

for file in *kraken.fa
do
kraken2-build --add-to-library $file --db Foxysporum
done

centrifuge-download -o library -d fungi -t taxon.oxysporum.txt refseq



srun --partition long --mem-per-cpu 10G --cpus-per-task 10 --pty bash
mkdir download
downloadRefSeq.pl --seqencesOutDirectory download/refseq --targetBranches fungi --taxonomyOutDirectory download/taxonomy




```bash
for Read1 in $(ls raw_data/S3_D6_readfish1/S3_D6_readfish_pass.fastq.gz); do
Out=$(echo $Read1 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $Out
Database=Foxysporum
OutDir=analysis/Kraken/$Out
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
sbatch $ProgDir/kraken_long_reads.sh $Read1 $Database $OutDir
done

for Read1 in $(ls raw_data/S3D6-2_pass.fastq.gz); do
Out=$(echo $Read1 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $Out
Database=Foxysporum
OutDir=analysis/Kraken/$Out
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
sbatch $ProgDir/kraken_long_reads.sh $Read1 $Database $OutDir
done

for Read1 in $(ls raw_data/S3D6_pass.fastq.gz); do
Out=$(echo $Read1 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $Out
Database=Foxysporum
OutDir=analysis/Kraken/$Out
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
sbatch $ProgDir/kraken_long_reads.sh $Read1 $Database $OutDir
done

for bar in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
Read1=raw_data/S3D6-BC11_BC12_No_Selection/"$bar"_pass.fastq.gz
Out=$(echo $Read1 | rev | cut -f1 -d '/' | rev | sed 's/_pass.fastq.gz//g')
echo $Out
Database=Foxysporum
OutDir=analysis/Kraken/S3D6-BC11_BC12_No_Selection/$Out
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Metagenomics
sbatch $ProgDir/kraken_long_reads.sh $Read1 $Database $OutDir
done



```

perl buildDB.pl --DB databases/FoxysporumDB --FASTAs database4kraken --taxonomy database4kraken/Foxysporum/taxonomy


## minimap2


```bash
    for Reference in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa); do
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