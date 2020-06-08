# Assembly qc tools

Tools used in the quality control and edition of genome assemblies

1. Quast. Quality assessment tool

2. Kat

3. 

4.  


## Quast

Produce genome assemblies statistics

### Requirements

```bash
# Python 2.7 required
conda install quast
```

### Typical run
```bash
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    for Assembly in $(ls path/to/genome/*.fa); do
        OutDir=$(dirname $Assembly)
        sbatch $ProgDir/quast.sh $Assembly $OutDir
    done
```