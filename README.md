# LDER
We design a method for heritability estimation, namely LD Eigenvalue Regression (LDER), which extends the LDSC method and provides more accurate estimates of heritability and confounding inflation.

Citation:

Shuang Song, Wei Jiang, Yiliang Zhang, Lin Hou, and Hongyu Zhao. Leveraging LD eigenvalue regression to improve the estimation of SNP heritability and confounding inflation. Submitted.

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
Users could either specify their own LD reference files with plink bfile format (.bim, .fam, .bed), or use the pre-computed LD information. We provide two examples here.

NOTE: We suggest users use plink bfile as the input, because the different numbers of SNPs in GWAS and in the reference panel may lead to a slight difference in the LD shrinkage.

### Example 1: Use plink bfile as the input (recommended)
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

### Example 2: Use the pre-computed LD information

The pre-computed LD information of 276,050 UK Biobank European individuals can be downloaded by

`wget -O LD.shrink.zip https://cloud.tsinghua.edu.cn/f/abf1020acb9c435eaa13/?dl=1 --no-check-certificate`

`wget -O LD.zip https://cloud.tsinghua.edu.cn/f/d93a4a7013fe461aa9fc/?dl=1 --no-check-certificate`

`unzip LD.shrink.zip`

`unzip LD.zip`

## Estimation of heritability and inflation factor
The main funcion can be run with:

```r
runLDER(assoc=GWAS_SUMMARY_STATISTICS (required), 
	n.gwas=SAMPLE_SIZE_OF_GWAS (required), 
	path=OUTPUT_DIR (required),
	LD.insample=IN_SAMPLE_LD (T/F, required),
	n.ld=SAMPLE_SIZE_OF_LD_REF (required), 
	ethnic=ETHNIC (optional),
	method=METHOD (default='lder')
	cores=NUMBER_OF_CORES (optional),
	a=INFLATION_FACTOR (optional))
```
- GWAS_SUMMARY_STATISTICS (required): GWAS summary statistics, need to include `snp`, `chr`, `a0`, `a1`, `z` (header is necessary)

- n.gwas (required): The sample size of the GWAS summary statistics

- OUTPUT_DIR (required): The output path (Note that the path should be SAME with that used in function `generateLD`)

- IN_SAMPLE_LD (required): T/F, whether the LD reference is estimated with the GWAS cohort (T) or external reference panel (e.g. 1000 Genome Project: F)

- SAMPLE_SIZE_OF_LD_REF (required): The sample size of the LD reference (e.g., 489 for 1000G)

- ETHNIC (optional): Ethnic of the GWAS cohort; 'eur' for European ancestry.

- METHOD (optional): Default='lder'. We also provide a choice of 'both', which outputs the results for both LDER and LDSC.

- NUMBER_OF_CORES (optional): The number of cores for computation in parallel.

- INFLATION_FACTOR (optional): Pre-specified inflation factor, default=NULL.



## Output

If `method='lder'`, the `runLDER` function returns a list with 4 elements:

`h2`: Estimated heritability by LDER

`inf`: Estimated inflation factor by LDER

`h2.se`: The standard error of estimated heritability estimated with block-jackknife.

`inf.se`: The standard error of estimated inflation factor estimated with block-jackknife.

If `method='both'`, the `runLDER` function returns a list containing the results of both LDER and LDSC.


## A Simplified Pipeline
Download a sample GWAS summary statistics:

$ wget -O gwas_sample.txt https://cloud.tsinghua.edu.cn/f/828ab71c87d84dd28d47/?dl=1

Download 1000G LD reference:

$ wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz

$ tar -xvzf 1000G_Phase3_plinkfiles.tgz


Run with R:

```r
devtools::install_github('shuangsong0110/LDER')
library(LDER)
library(data.table)
path0 <- getwd()
assoc <- fread('gwas_sample.txt')
for(chr in 1:22){
    generateLD(assoc, path = path0, bfile_path = paste0('/1000G_EUR_Phase3_plink/1000G.EUR.QC.', chr))
}
res <- runLDER(assoc, n.gwas=2e4, path=path0, LD.insample=F, ethnic='eur', n.ld=489, cores=10, method='lder', a=NULL)

```


## Maintainer

Please contact Shuang Song (song-s19@mails.tsinghua.edu.cn) if there are any problems or questions.


