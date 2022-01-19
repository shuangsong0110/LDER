#' @title Calculate LD files
#' @description Run LDER
#' @param assoc GWAS summary statistics, need to include snp, chr, a0, a1, z (header is necessary)
#' @param path The output path
#' @param bfile_path The path to the plink bfiles for LD estimation
#' @param cores The number of cores for computation in parallel
#' @param ethnic Ethnic of the GWAS cohort; 'eur' for European ancestry
#' @param plink_path The path to the plink software. If not specified, the function will use the default path (which plink)
#' @param python_path The path to the python software. If not specified, the function will use the default path (which python)
#' @import  data.table stats utils
#' @export


generateLD <- function(assoc, path, bfile_path, cores=10, ethnic='eur', plink_path=NULL, python_path=NULL){
  dir.create(path,recursive = T)
  setwd(path)
  bim <- fread(paste0(bfile_path,'.bim'))
  if(length(intersect(bim$V2,assoc$snp))>0){
  if(is.null(python_path)){
    temp <- system('which python')
    if(temp==0){
      python <- system('which python',intern=TRUE)
    }else{
      warning('Please install python or provide the exact path to python software')
    }
  }else{
    python <- python_path
  }
  if(is.null(plink_path)){
    temp <- system('which plink')
    if(temp==0){
      plink <- system('which plink',intern=TRUE)
    }else{
      warning('Please install plink or provide the exact path to plink software')
    }
  }else{
    plink <- plink_path
  }
  if(!file.exists('./plinkLD/plinkLD.py')){
    system('wget -O plinkLD.zip https://cloud.tsinghua.edu.cn/f/3f96074ee7ee436895ac/?dl=1 --no-check-certificate')
  system('unzip plinkLD.zip')
  system('rm plinkLD.zip')
  }
  write.table(assoc$snp,'snp.txt',quote=F,row.names=F,col.names=F)
  system(paste0(plink,' --bfile ',bfile_path,' --extract snp.txt --make-bed --out geno --noweb'))
  system(paste0(python,' ./plinkLD/plinkLD.py --bfile geno  --block ./plinkLD/ldetect-data/fourier_ls-all.bed --snplist ./plinkLD/w_hm3.snplist --output ',
  path,'/LD.shrink/ --thread ',cores,' --method LW'))
  }
}
