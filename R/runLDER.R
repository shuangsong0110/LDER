#' @title Main function
#' @description Run LDER
#' @param assoc GWAS summary statistics, need to include snp, chr, a0, a1, z (header is necessary)
#' @param n.gwas The sample size of the GWAS summary statistics
#' @param path The output path
#' @param LD.insample T/F, whether the LD reference is estimated with the GWAS cohort (T) or external reference panel (e.g. 1000 Genome Project) (F)
#' @param ethnic Ethnic of the GWAS cohort; 'eur' for European ancestry
#' @param n.ld The sample size of the LD reference
#' @param cores The number of cores for computation in parallel
#' @param method 'lder', 'ldsc', or 'both'
#' @param a Pre-specified inflation factor, default=NULL
#' @import  data.table stats utils
#' @export
#'
#'
runLDER <- function(assoc, n.gwas, path, LD.insample=F,  n.ld, ethnic='eur',method='lder', cores=10, a=NULL, rough=T, twostage=T, type='jack', n.bs=1000){
  dir.create(path,recursive = T)

  setwd(path)
  rough <- !LD.insample
  bed <- fread('./plinkLD/ldetect-data/fourier_ls-all.bed')
  library(parallel)
  ldpath <- ldpath.shrink <- paste0(path,'/LD.shrink/')
  stats <- mclapply(0:(nrow(bed)-1),get.stats,assoc=assoc,ldpath=ldpath,ldpath.shrink=ldpath.shrink,mc.cores=cores,n.ld=n.ld)
  if(method=='lder'){
    res <- lder(stats=stats,n.gwas=n.gwas,a=a,rough=rough,cores=cores,twostage=twostage,type=type,n.bs=n.bs)
    return(res)
  }else if(method=='ldsc'){
    res <- ldsc(stats=stats,n.gwas=n.gwas,a=a,cores=cores,twostage=twostage,type=type,n.bs=n.bs)
    return(res)
  }else if(method=='both'){
    res1 <- lder(stats=stats,n.gwas=n.gwas,a=a,rough=rough,cores=cores,twostage=twostage,type=type,n.bs=n.bs)
    res2 <- ldsc(stats=stats,n.gwas=n.gwas,a=a,cores=cores,twostage=twostage,type=type,n.bs=n.bs)
    return(list(lder=res1,ldsc=res2))
  }
}



