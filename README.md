# LDER
We design a method for heritability estimation, namely LD Eigenvalue Regression (LDER), which extends the LDSC method and provides more accurate estimates of heritability and confounding inflation.

## Table of contents
* [Install](#install)
* [LD prepared](#ld-prepared)
* [Estimation of heritability and inflation factor](#estimation-of-heritability-and-inflation-factor)
* [Output](#output)
* [A Simplified Pipeline](#a-simplified-pipeline)

## Install
LDER is an R package which can be installed using the command:
```r
devtools::install_github('shuangsong0110/LDER')
```

## LD prepared
We provide a function `plinkLD.py` for efficient LD information extraction and shrinkage based on Python. 
Users could specify their own LD reference files with plink bfile format (.bim, .fam, .bed).
The 1000 Genome Project reference panel (hg19) can be downloaded by:

`wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz`

`tar -xvzf 1000G_Phase3_plinkfiles.tgz`

```r
library(LDER)
generateLD(assoc=GWAS_SUMMARY_STATISTICS (required), 
          path=OUTPUT_DIR (required),
          bfile_path=PATH_TO_LD_REFERENCE (required),
          cores=NUMBER_OF_CORES (optional),
          ethnic=ETHNIC (optional),
          plink_path=PATH_TO_PLINK_SOFTWARE (optional),
          python_path=PATH_TO_PYTHON_SOFTWARE (optional))                    
```
- GWAS_SUMMARY_STATISTICS (required): GWAS summary statistics, need to include `snp`, `chr`, `a0`, `a1`, `z` (header names are necessary)

- OUTPUT_DIR (required): The output path

- PATH_TO_LD_REFERENCE (required): The LD reference plink bfile. If the files are divided into chromosomes, the function should be run for each of the file.

- NUMBER_OF_CORES (optional): The number of cores for computation in parallel.

- ETHNIC (optional): Ethnic of the GWAS cohort; 'eur' for European ancestry.

- PATH_TO_PLINK_SOFTWARE (optional): The path to the plink software. If not specified, the function will use the default path (system("which plink"))

- PATH_TO_PYTHON_SOFTWARE (optional): The path to the python software. If not specified, the function will use the default path (system("which python"))








