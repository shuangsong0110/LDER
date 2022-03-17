### input: z score, a0,a1, LD files (plinkLD output)
source('/home/songs/genomic_control/ukb/20210204/ldsc_new.R')
library(data.table)

agtc <- function(a1,a2,b1,b2){
  sig <- rep(1,length(a1))
  for(i in 1:length(a1)){
    if((a1[i]=='A')&(a2[i]=='T')){
      sig[i] <- 0
    }else if((a1[i]=='T')&(a2[i]=='A')){
      sig[i] <- 0
    }else if((a1[i]=='C')&(a2[i]=='G')){
      sig[i] <- 0
    }else if((a1[i]=='G')&(a2[i]=='C')){
      sig[i] <- 0
    }else if((a1[i]==b1[i])&(a2[i]==b2[i])){
      sig[i] <- 1
    }else if((a1[i]==b2[i])&(a2[i]==b1[i])){
      sig[i] <- -1
    }else{
      if(b1[i]=="A"){temp1 <- "T"}
      if(b1[i]=="T"){temp1 <- "A"}
      if(b1[i]=="G"){temp1 <- "C"}
      if(b1[i]=="C"){temp1 <- "G"}
      if(b2[i]=="A"){temp2 <- "T"}
      if(b2[i]=="T"){temp2 <- "A"}
      if(b2[i]=="G"){temp2 <- "C"}
      if(b2[i]=="C"){temp2 <- "G"}
      if((a1[i]==temp1)&(a2[i]==temp2)){
        sig[i] <- 1
      }else if((a1[i]==temp2)&(a2[i]==temp1)){
        sig[i] <- -1
      }else{ sig[i] <- 0}
    }
  }
  kk <- length(which(sig==0))
  #print(paste0('Remove ',kk,' SNPs due to ambiguous alleles'))
  return(sig)
}
bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm=T)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}
# assoc:: z, chr, snp, a0, a1

get.stats <- function(j,assoc,ldpath,ldpath.shrink,n.ld=1e4){
  ## read info
  if(is.null(ldpath.shrink)){
    ldpath.shrink <- ldpath
  }
  if(file.exists(paste0(ldpath.shrink, '/SNPINFO', j, '.txt'))){
    ld.info <- fread(paste0(ldpath.shrink, '/SNPINFO', j, '.txt'))
    ld.info$order <- 1:nrow(ld.info)
    if(length(which(is.na(assoc$z)))>0){
      assoc <- assoc[-which(is.na(assoc$z)),]
    }
    assoc.sub <- assoc[which(assoc$chr==ld.info$CHR[1]),]
    assoc.new <- merge(assoc.sub,ld.info,by.x='snp',by.y='SNP')
    if(dim(assoc.new)[1]>0){
      sign <- agtc(assoc.new$a0, assoc.new$a1, assoc.new$A1, assoc.new$A2)
      assoc.new$z <- assoc.new$z*sign
      assoc.new <- assoc.new[which(sign!=0),]
      if(dim(assoc.new)[1]>0){
        #print(dim(assoc.new)[1])
        assoc.new <- assoc.new[order(assoc.new$BP),]
        ## read LD shrink files 
        if(nrow(ld.info)==1){
          ld.shrink <- 1
        }else{
          #if(file.exists(paste0(ldpath.shrink, '/LD', j, '.txt'))){
          ld.shrink <- as.matrix(fread(paste0(ldpath.shrink, '/LD', j, '.txt')))
          R0 <- ld.shrink[assoc.new$order,assoc.new$order]
          R0[which(is.na(R0))] <- 0
        }
        temp <- eigen(R0,symmetric = T)
        #eff.num <- length(which(temp$values>0))
        temp$values[temp$values<1e-6] <- 0
        U <- temp$vectors
        V <- diag(temp$values)
        V.inv <-  1/V
        V.inv[which(V.inv==Inf)] <- 0
        #eff.num <- length(which(V>0))
        eff.num <- dim(R0)[1]
        if(is.null(eff.num)){eff.num <- 1}
        eigen.mat<- sqrt(V.inv[1:eff.num,1:eff.num])%*%(t(U)[1:eff.num,])
        lam <- diag(V)[1:eff.num] ## export
        z <- assoc.new$z
        x <- eigen.mat%*%z ## export
        if(nrow(ld.info)==1){
          R0 <- 1
        }else{
          ld <- as.matrix(fread(paste0(ldpath, '/LD', j, '.txt')))
          R0 <- ld[assoc.new$order,assoc.new$order]
          R0[which(is.na(R0))] <- 0
        }
        ldsc <- ldscore(R0,N=n.ld) ## export
        return(list(x=x,z=z,lam=lam,ldsc=ldsc))
      }else{
        print(paste0('0 SNPs in LD ',j))
      }
    }else{
      print(paste0('0 SNPs in LD ',j))
    } 
  }else{
    print(paste0(ldpath.shrink, '/SNPINFO', j, '.txt'))
    print(paste0('SNPINFO',j, ' does not exist'))
  }
  #}
}
#get.stats(188,assoc,ldpath,ldpath.shrink)
get.res <- function(x1,lam1,n.gwas,a,rough,twostage){
  m <- length(x1)
  if(twostage){
    newGC <-  1+n.gwas*calH2.new1(x1, lam1,   N=n.gwas,a=NULL, rough=F)$a
    #print(newGC)
    if(newGC<1.05){
      s2thld <- 3.5*sqrt(m/n.gwas)+2.5 
    }else{
      s2thld <- quantile(x1^2, 1-5/m)
    }
    idx <- which(x1^2<=s2thld)
    temp1 <-  calH2.new1(x1[idx], lam1[idx],   N=n.gwas,a=a, rough=rough)$a
    a <- max(temp1,0)
    print(a*n.gwas)
    h2 <- calH2.new1(x1, lam1,   N=n.gwas, a=max(temp1,0), rough=rough)$h2
  }else{
      a <- calH2.new1(x1, lam1,   N=n.gwas, a=a, rough=rough)$a
      a <- max(a,0)
      h2 <- calH2.new1(x1, lam1,   N=n.gwas, a=a, rough=rough)$h2
      
  }
  return(list(a=a,h2=h2))
}

get.res.ldsc <- function(z1,ldsc1,n.gwas,a,twostage){
  m <- length(z1)
  s1thld <- 30
  if(twostage){
    a <-  calH2(z1[(z1^2)<s1thld], ldsc1[(z1^2)<s1thld],   N=n.gwas,a=NULL)$a
    h2 <- calH2(z1, ldsc1,   N=n.gwas,a=a)$h2
  }else{
    h2 <-  calH2(z1, ldsc1,   N=n.gwas,a=a)$h2
    a <- calH2(z1, ldsc1,   N=n.gwas,a=a)$a
  }
  return(list(a=a,h2=h2))
}

cvt.h2 <- function(K,P,h2){
  # K: prevalence
  # P: proportion of cases
  zv <- dnorm(qnorm(K))
  h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
  return(h2_liab)
}

lder <- function(stats,n.gwas,a=NULL,rough=F,cores=20,twostage=T,type='jack',n.bs=1000){
  library(parallel)
  #stats <- mclapply(0:1702,get.stats,assoc=assoc,ldpath=ldpath,ldpath.shrink=ldpath.shrink,mc.cores=cores)
  temp <- unlist(lapply(stats,length))
  stats[temp==1] <- NULL
  x1 <- unlist(sapply(stats,'[[',1))
  lam1 <- unlist(sapply(stats,'[[',3))
  do.jack <- function(stats,j,cores,n.gwas,a,rough,twostage){
    stats1 <- stats[-j]
    x1 <- unlist(sapply(stats1,'[[',1))
    lam1 <- unlist(sapply(stats1,'[[',3))
    
    res <- get.res(x1,lam1,n.gwas=n.gwas,a=a,rough=rough,twostage=twostage)
    res$h2 <- res$h2/length(x1)
    return(res)
  }
  do.boot <- function(stats,j,cores,n.gwas,a,rough,twostage){
    #print(j)
    n.block <- length(stats)
    set.seed(j)
    ind.bs <- sample(n.block,n.block,replace = T)
    stats1 <- stats[ind.bs]
    x1 <- unlist(sapply(stats1,'[[',1))
    z1 <- unlist(sapply(stats1,'[[',2))
    lam1 <- unlist(sapply(stats1,'[[',3))
    res <- get.res(x1,lam1,n.gwas=n.gwas,a=a,rough=rough,twostage=twostage)
    res$h2 <- res$h2/length(x1)
    
    return(res)
  }
  
  res <- get.res(x1,lam1,n.gwas=n.gwas,a=a,rough=rough,twostage=twostage)
  if(type=='jack'){
    ### jackknife
    temp.jack <- mclapply(1:length(stats),do.jack ,stats=stats,
                          n.gwas=n.gwas,a=a,rough=rough,twostage=twostage,mc.cores =cores,mc.preschedule=F)
    hh2 <-  unlist(sapply(temp.jack,'[[','h2'))
    inff2 <- unlist(sapply(temp.jack,'[[','a'))
    n.block <- length(stats)
    
    h33.jack <- sqrt((n.block-1)/n.block*sum((hh2-mean(hh2))^2))*length(x1)
    a33.jack <- sqrt((n.block-1)/n.block*sum((inff2-mean(inff2))^2))
    
    return(list(h2=res$h2,inf=res$a*n.gwas+1,h2.sd=h33.jack,inf.sd=a33.jack*n.gwas))
  }else if(type=='boot'){
    temp.jack <- mclapply(1:n.bs,do.boot ,stats=stats,
                          n.gwas=n.gwas,a=a,rough=rough,twostage=twostage,mc.cores =cores,mc.preschedule=F)
    hh2 <-  unlist(sapply(temp.jack,'[[','h2'))
    inff2 <- unlist(sapply(temp.jack,'[[','a'))
    n.block <- length(stats)
    h33.jack <- sqrt(1/(n.bs)*sum((hh2-mean(hh2))^2))*length(x1)
    a33.jack <- sqrt(1/(n.bs)*sum((inff2-mean(inff2))^2))
    return(list(h2=res$h2,inf=res$a*n.gwas+1,h2.sd=h33.jack,inf.sd=a33.jack*n.gwas))
  }else if(type=='both'){
    temp.jack <- mclapply(1:length(stats),do.jack ,stats=stats,
                          n.gwas=n.gwas,a=a,rough=rough,twostage=twostage,mc.cores =cores,mc.preschedule=F)
    hh2 <-  unlist(sapply(temp.jack,'[[','h2'))
    inff2 <- unlist(sapply(temp.jack,'[[','a'))
    n.block <- length(stats)
    
    h33.jack <- sqrt((n.block-1)/n.block*sum((hh2-mean(hh2))^2))*length(x1)
    a33.jack <- sqrt((n.block-1)/n.block*sum((inff2-mean(inff2))^2))
    
    temp.jack <- mclapply(1:n.bs,do.boot ,stats=stats,
                          n.gwas=n.gwas,a=a,rough=rough,twostage=twostage,mc.cores =cores,mc.preschedule=F)
    hh2 <-  unlist(sapply(temp.jack,'[[','h2'))
    inff2 <- unlist(sapply(temp.jack,'[[','a'))
    n.block <- length(stats)
    h33.boot <- sqrt(1/(n.bs)*sum((hh2-mean(hh2))^2))*length(x1)
    a33.boot <- sqrt(1/(n.bs)*sum((inff2-mean(inff2))^2))
    return(list(h2=res$h2,inf=res$a*n.gwas+1,h2.sd.jack=h33.jack,inf.sd.jack=a33.jack*n.gwas,
                h2.sd.boot=h33.boot,inf.sd.boot=a33.boot*n.gwas ))
  }else if(type=='none'){
    return(list(h2=res$h2,inf=res$a*n.gwas+1 ))
  }
}

ldsc <- function(stats,n.gwas,a=NULL,cores=20,twostage=T,type='jack',n.bs=1000){
  library(parallel)
  #stats <- mclapply(0:1702,get.stats,assoc=assoc,ldpath=ldpath,ldpath.shrink=ldpath.shrink,mc.cores=cores)
  temp <- unlist(lapply(stats,length))
  stats[temp==1] <- NULL
  z1 <- unlist(lapply(stats,'[[','z'))
  ldsc1 <- unlist(sapply(stats,'[[','ldsc'))
  do.jack.ldsc <- function(stats,j,cores,n.gwas,a,rough,twostage){
    stats1 <- stats[-j]
    z1 <- unlist(sapply(stats1,'[[','z'))
    ldsc1 <- unlist(sapply(stats1,'[[','ldsc'))
    res <- get.res.ldsc(z1,ldsc1,n.gwas=n.gwas,a=a,twostage=twostage)
    res$h2 <- res$h2/length(z1)
    return(res)
  }
  do.boot.ldsc <- function(stats,j,cores,n.gwas,a,rough,twostage){
    #print(j)
    n.block <- length(stats)
    set.seed(j)
    ind.bs <- sample(n.block,n.block,replace = T)
    stats1 <- stats[ind.bs]
    z1 <- unlist(sapply(stats1,'[[','z'))
    ldsc1 <- unlist(sapply(stats1,'[[','ldsc'))
    res <- get.res.ldsc(z1,ldsc1,n.gwas=n.gwas,a=a,twostage=twostage)
    res$h2 <- res$h2/length(z1)
    return(res)
  }
  
  res <- get.res.ldsc(z1,ldsc1,n.gwas=n.gwas,a=a,twostage=twostage)
  if(type=='jack'){
    ### jackknife
    temp.jack <- mclapply(1:length(stats),do.jack.ldsc ,stats=stats,
                          n.gwas=n.gwas,a=a,rough=rough,twostage=twostage,mc.cores =cores,mc.preschedule=F)
    hh2 <-  unlist(sapply(temp.jack,'[[','h2'))
    inff2 <- unlist(sapply(temp.jack,'[[','a'))
    n.block <- length(stats)
    
    h33.jack <- sqrt((n.block-1)/n.block*sum((hh2-mean(hh2))^2))*length(z1)
    a33.jack <- sqrt((n.block-1)/n.block*sum((inff2-mean(inff2))^2))
    
    return(list(h2=res$h2,inf=res$a*n.gwas+1,h2.sd=h33.jack,inf.sd=a33.jack*n.gwas))
  }else if(type=='boot'){
    temp.jack <- mclapply(1:n.bs,do.boot.ldsc ,stats=stats,
                          n.gwas=n.gwas,a=a,rough=rough,twostage=twostage,mc.cores =cores,mc.preschedule=F)
    hh2 <-  unlist(sapply(temp.jack,'[[','h2'))
    inff2 <- unlist(sapply(temp.jack,'[[','a'))
    n.block <- length(stats)
    h33.jack <- sqrt(1/(n.bs)*sum((hh2-mean(hh2))^2))*length(z1)
    a33.jack <- sqrt(1/(n.bs)*sum((inff2-mean(inff2))^2))
    return(list(h2=res$h2,inf=res$a*n.gwas+1,h2.sd=h33.jack,inf.sd=a33.jack*n.gwas))
    
  }else if (type=='both'){
    temp.jack <- mclapply(1:length(stats),do.jack.ldsc ,stats=stats,
                          n.gwas=n.gwas,a=a,rough=rough,twostage=twostage,mc.cores =cores,mc.preschedule=F)
    hh2 <-  unlist(sapply(temp.jack,'[[','h2'))
    inff2 <- unlist(sapply(temp.jack,'[[','a'))
    n.block <- length(stats)
    
    h33.jack <- sqrt((n.block-1)/n.block*sum((hh2-mean(hh2))^2))*length(z1)
    a33.jack <- sqrt((n.block-1)/n.block*sum((inff2-mean(inff2))^2))
    
    temp.jack <- mclapply(1:n.bs,do.boot.ldsc ,stats=stats,
                          n.gwas=n.gwas,a=a,rough=rough,twostage=twostage,mc.cores =cores,mc.preschedule=F)
    hh2 <-  unlist(sapply(temp.jack,'[[','h2'))
    inff2 <- unlist(sapply(temp.jack,'[[','a'))
    n.block <- length(stats)
    h33.boot <- sqrt(1/(n.bs)*sum((hh2-mean(hh2))^2))*length(z1)
    a33.boot <- sqrt(1/(n.bs)*sum((inff2-mean(inff2))^2))
    return(list(h2=res$h2,inf=res$a*n.gwas+1,h2.sd.jack=h33.jack,inf.sd.jack=a33.jack*n.gwas,
                h2.sd.boot=h33.boot,inf.sd.boot=a33.boot*n.gwas ))
  }else if(type=='none'){
    return(list(h2=res$h2,inf=res$a*n.gwas+1))
  }
}

ld.all <- function(assoc,n.gwas,ldpath,ldpath.shrink=NULL,method='lder',
                   a=NULL,rough=F,cores=20,twostage=T,type='jack',n.bs=1000,n.ld=276050){
  library(parallel)
  stats <- mclapply(0:1702,get.stats,assoc=assoc,ldpath=ldpath,ldpath.shrink=ldpath.shrink,mc.cores=cores,n.ld=n.ld)
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



