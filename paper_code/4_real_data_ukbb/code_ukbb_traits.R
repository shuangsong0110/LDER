
#######
path <- '/home/songs/UKB_summary/20210211/'
setwd(path)
source('/home/songs/genomic_control/20210714_code.R')

library(data.table)
dat <- fread('0224dat.analyse.txt',sep='\t',head=T)
## File: file name
## 
dat <- as.data.frame(dat)


#get.stats(81,assoc,ldpath,ldpath.shrink, path.pheno, n.ld)
cores <- 45
n.ld <- 276050

run.res <- function(w,dat){
  setwd(path)
  pheno <- dat$pheno.code[w]
  print(pheno)
  path.res <- paste0(path,'/res/',pheno)
  fn2 <- paste0('./clean/',pheno,'.summary.txt')
  if((!file.exists(paste0('./res_ukb/newinf_quantitative/',pheno,'_res_220212.RData')))&
     (!file.exists(paste0('./res_ukb/newinf_qualitative/',pheno,'_res_220212.RData')))&
     (file.exists(fn2))&
     (file.info(fn2)$size>0)){
    assoc0 <- as.data.frame(fread(fn2))
    assoc <- assoc0[,c('V4','V2','a0','a1','chr','z','N')]
    colnames(assoc)[1:2] <- c('bp','snp')
    
    library(data.table)
    ldpath <- '/home/songs/genomic_control/ukb/20210204/ukb/LD/'
    ldpath.shrink <-  '/home/songs/genomic_control/ukb/20210204/ukb/LD.shrink/'
    res <- ld.all(assoc,assoc$N[1],ldpath=ldpath,ldpath.shrink=ldpath.shrink,method='both',
                  type='jack',rough=F,twostage=T,cores=10)
     if(dat$variable_type[w]!='binary'){
      setwd('/home/songs/UKB_summary/20210211/res_ukb/newinf_quantitative')
      save(res,file=paste0(pheno,'_res_220212.RData'))
      print('done')
    }else{
      setwd('/home/songs/UKB_summary/20210211/res_ukb/newinf_qualitative')
      ncase <- dat$n_cases[w]
      ncontrol <- dat$n_controls[w]
      K <- P <- ncase/(ncase+ncontrol)
      res$lder$h2.lia <- cvt.h2(K,P,res$lder$h2)
      res$lder$h2.sd.lia <- cvt.h2(K,P,res$lder$h2.sd)
      res$ldsc$h2.lia <- cvt.h2(K,P,res$ldsc$h2)
      res$ldsc$h2.sd.lia <- cvt.h2(K,P,res$ldsc$h2.sd)
      print(res)
      save(res,file=paste0(pheno,'_res_220212.RData'))
      print('done')
    }
  }
  
}
run.res(2,dat)
library(parallel)
quan.trait <- which(dat$variable_type!='binary')
qual.trait <- which(dat$variable_type=='binary')
res <- mclapply(quan.trait,run.res,dat=dat,mc.cores=5)
res <- mclapply(qual.trait,run.res,dat=dat,mc.cores=5)


res <- mclapply(1:nrow(dat),run.res,dat=dat,mc.cores=30)


for(w in 1:nrow(dat)){
  setwd(path)
  print(w)
  pheno <- dat$pheno.code[w]
  path.res <- paste0(path,'/res/',pheno)
  fn2 <- paste0('./clean/',pheno,'.summary.txt')
  if((!file.exists(paste0('./res_ukb/quantitative/',pheno,'_res.RData')))&
     (!file.exists(paste0('./res_ukb/qualitative/',pheno,'_res.RData')))&
     (file.exists(fn2))&
     (file.info(fn2)$size>0)){
    assoc0 <- as.data.frame(fread(fn2))
    assoc <- assoc0[,c('V4','V2','a0','a1','chr','z','N')]
    colnames(assoc)[1:2] <- c('bp','snp')
    
    library(data.table)
    ldpath <- '/home/songs/genomic_control/ukb/20210204/ukb/LD/'
    ldpath.shrink <-  '/home/songs/genomic_control/ukb/20210204/ukb/LD.shrink/'
    res <- ld.all(assoc,assoc$N[1],ldpath=ldpath,ldpath.shrink=ldpath.shrink,method='both',
                  type='jack',rough=F,twostage=T,cores=40)
    
    if(dat$variable_type[w]!='binary'){
      setwd('/home/songs/UKB_summary/20210211/res_ukb/quantitative')
      save(res,file=paste0(pheno,'_res.RData'))
    }else{
      setwd('/home/songs/UKB_summary/20210211/res_ukb/qualitative')
      ncase <- dat$n_cases[w]
      ncontrol <- dat$n_controls[w]
      K <- P <- ncase/(ncase+ncontrol)
      res$lder$h2.lia <- cvt.h2(K,P,res$lder$h2)
      res$lder$res <- cvt.h2(K,P,res$lder$h2.sd)
      res$ldsc$h2.lia <- cvt.h2(K,P,res$ldsc$h2)
      res$ldsc$h2.sd.lia <- cvt.h2(K,P,res$ldsc$h2.sd)
      print(res)
      save(res,file=paste0(pheno,'_res.RData'))
    }
  }
}
if(T){
  ## save prevalence
  dat <- as.data.frame(dat)
  dat <- dat[!duplicated(dat$pheno.code),]
  
  prev <- c()
  jj <- 1
  for(w in 1:nrow(dat)){
    print(jj)
    setwd(path)
    print(w)
    pheno <- dat$pheno.code[w]
    path.res <- paste0(path,'/res/',pheno)
    fn2 <- paste0('./clean/',pheno,'.summary.txt')
    if((file.exists(fn2))&
       (file.info(fn2)$size>0)){
      library(data.table)
      if(dat$variable_type[w]!='binary'){
        setwd('/home/songs/UKB_summary/20210211/res_ukb/quantitative')
        prev[jj] <- 1
      }else{
        setwd('/home/songs/UKB_summary/20210211/res_ukb/qualitative')
        ncase <- dat$n_cases[w]
        ncontrol <- dat$n_controls[w]
        K <- P <- ncase/(ncase+ncontrol)
        prev[jj] <- K
      }
      jj <- jj+1
    }
    
  }
  setwd("/home/songs/UKB_summary/20210211/res_ukb/")
  save(prev,file='prev.RData')
}


##### save result

#######
path <- '/home/songs/UKB_summary/20210211/'
setwd(path)

library(data.table)
dat <- fread('0224dat.analyse.txt',sep='\t',head=T)
## File: file name
## 
dat <- as.data.frame(dat)
dat <- dat[!duplicated(dat$pheno.code),]
dat.quan <- dat[dat$variable_type!='binary',]
dat.qual <- dat[dat$variable_type=='binary',]
if(T){
  dat <- dat.quan
  h2.ldsc <- h2.lder <- h2.ldsc.sd <- h2.lder.sd <- inf.lder <- 
    inf.ldsc <- inf.lder.sd <- inf.ldsc.sd <- c()
  for(w in 1:nrow(dat)){
    setwd(path)
    print(w)
    pheno <- dat$pheno.code[w]
    setwd('/home/songs/UKB_summary/20210211/res_ukb/newinf_quantitative')
    load(paste0(pheno,'_res_220212.RData'))
    h2.lder <- c(h2.lder,res$lder$h2)
    h2.lder.sd <- c(h2.lder.sd,res$lder$h2.sd)
    inf.lder <- c(inf.lder,res$lder$inf)
    inf.lder.sd <- c(inf.lder.sd,res$lder$inf.sd)
    
    h2.ldsc <- c(h2.ldsc,res$ldsc$h2)
    h2.ldsc.sd <- c(h2.ldsc.sd,res$ldsc$h2.sd)
    inf.ldsc <- c(inf.ldsc,res$ldsc$inf)
    inf.ldsc.sd <- c(inf.ldsc.sd,res$ldsc$inf.sd)
    
  }
  res <- data.frame(Trait=dat$pheno.code,
                    h2.LDER= paste0(format(round(h2.lder,3), nsmall = 3),' (',format(round(h2.lder.sd,3), nsmall = 3),')'),
                    h2.LDSC= paste0(format(round(h2.ldsc,3), nsmall = 3),' (',format(round(h2.ldsc.sd,3), nsmall = 3),')'),
                    inf.LDER= paste0(format(round(inf.lder,3), nsmall = 3),' (',format(round(inf.lder.sd,3), nsmall = 3),')'),
                    inf.LDSC= paste0(format(round(inf.ldsc,3), nsmall = 3),' (',format(round(inf.ldsc.sd,3), nsmall = 3),')'),
                    Description=dat$description)
  # res <- res[!duplicated(res$Trait),]
  setwd('/home/songs/UKB_summary/20210211/res_ukb')
  write.csv(res,'newinf_res_quan.csv',quote=T,row.names = F)
  #write.table(res,file="res_quan_new.csv",quote=F,row.names=F,col.names = T)
}



if(T){
  dat <- dat.qual
  h2.ldsc <- h2.lder <- h2.ldsc.sd <- h2.lder.sd <- inf.lder <- 
    inf.ldsc <- inf.lder.sd <- inf.ldsc.sd <- c()
  for(w in 1:nrow(dat)){
    setwd(path)
    print(w)
    pheno <- dat$pheno.code[w]
    setwd('/home/songs/UKB_summary/20210211/res_ukb/qualitative')
    load(paste0(pheno,'_res.RData'))
    h2.lder <- c(h2.lder,res$lder$h2.lia)
    h2.lder.sd <- c(h2.lder.sd,res$lder$h2.sd.lia)
    inf.lder <- c(inf.lder,res$lder$inf)
    inf.lder.sd <- c(inf.lder.sd,res$lder$inf.sd)
    
    h2.ldsc <- c(h2.ldsc,res$ldsc$h2.lia)
    h2.ldsc.sd <- c(h2.ldsc.sd,res$ldsc$h2.sd.lia)
    inf.ldsc <- c(inf.ldsc,res$ldsc$inf)
    inf.ldsc.sd <- c(inf.ldsc.sd,res$ldsc$inf.sd)
    
  }
  res <- data.frame(Trait=dat$pheno.code,
                    h2.LDER= paste0(format(round(h2.lder,3), nsmall = 3),' (',format(round(h2.lder.sd,3), nsmall = 3),')'),
                    h2.LDSC= paste0(format(round(h2.ldsc,3), nsmall = 3),' (',format(round(h2.ldsc.sd,3), nsmall = 3),')'),
                    inf.LDER= paste0(format(round(inf.lder,3), nsmall = 3),' (',format(round(inf.lder.sd,3), nsmall = 3),')'),
                    inf.LDSC= paste0(format(round(inf.ldsc,3), nsmall = 3),' (',format(round(inf.ldsc.sd,3), nsmall = 3),')'),
                    Description=dat$description)
  setwd('/home/songs/UKB_summary/20210211/res_ukb')
  write.csv(res,'res_qual.csv',quote=T,row.names = F)
}






for(w in 1:nrow(dat)){
  setwd(path)
  print(w)
  pheno <- dat$pheno.code[w]
  N <- dat$n_non_missing[w]
  path.boot <- paste0(path,'/res.boot/',pheno)
  if(!file.exists(paste0(path.boot,'/res_1kgref.all.RData'))){
    res <- onekg.bs(path=path,pheno,n.gwas=N,n.bs=1000,cores=30)
    dir.create(path.boot,recursive=T)
    save(res,file=paste0(path.boot,'/res_1kgref.all.RData'))
  }
}




library(data.table)
dat=fread('result_type.csv')
length(which((dat$pv.lder<0.05/815)&(dat$pv.ldsc>=0.05/815)))

length(which((dat$pv.lder>=0.05/815)&(dat$pv.ldsc<0.05/815)))

