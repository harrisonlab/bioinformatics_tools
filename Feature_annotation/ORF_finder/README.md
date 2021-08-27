# Commands for the analysis of Turnip mosaic virus

### Extract seq from Whole_genomes

```bash
for genome in $(ls Whole_genomes/*seq); do
Strain=$(echo $genome | rev | cut -f1 -d '/' | rev | sed 's/.seq//g')
echo $Strain
echo ">"$Strain"_genome" > "$Strain"_temp.txt
OutDir=Genomes_edited/$Strain
mkdir -p $OutDir
less $genome | tail -2 > "$Strain"_temp.fa
cat "$Strain"_temp.txt "$Strain"_temp.fa > $OutDir/"$Strain".fa
rm "$Strain"_temp.fa
rm "$Strain"_temp.txt
done
```