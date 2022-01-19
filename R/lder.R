
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
    return(list(h2=res$h2,inf=res$a*n.gwas+1,h2.se=h33.jack,inf.se=a33.jack*n.gwas))
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

    return(list(h2=res$h2,inf=res$a*n.gwas+1,h2.se=h33.jack,inf.se=a33.jack*n.gwas))
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


