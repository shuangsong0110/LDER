if(T){
  source('/home/songs/genomic_control/20210714_code.R')
  
  for(trait in c('HEIGHT','BMI','HDL','LDL','TG','TC')){
    setwd('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/')
    dir.create(trait)
    print(trait)
    setwd(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait))
    load('lder.assoc.RData')
    ldpath <- '/home/songs/genomic_control/ukb/20210204/ukb/LD/'
    ldpath.shrink <- '/home/songs/genomic_control/ukb/20210204/ukb/LD.shrink/'
    
    res <- ld.all(assoc,assoc$N[1],ldpath=ldpath,ldpath.shrink=ldpath.shrink,method='lder',
                  type='jack',rough=F,twostage=T,cores=5)
    print(res)
    save(res,file='res_new_lder.RData')
  }
}

if(T){
  for(trait in c('T2D','ASTHMA','CAD','SCZ')){
    setwd('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/')
    dir.create(trait)
    print(trait)
    setwd(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait))
    load('lder.assoc.RData')
    ldpath <- '/home/songs/genomic_control/ukb/20210204/ukb/LD/'
    ldpath.shrink <- '/home/songs/genomic_control/ukb/20210204/ukb/LD.shrink/'
    
    res <- ld.all(assoc,assoc$N[1],ldpath=ldpath,ldpath.shrink=ldpath.shrink,method='lder',
                  type='jack',rough=F,twostage=T,cores=10)
    
    phen <- fread('plink.log',fill=T)
    ncase <- as.numeric(phen$V4[25])
    ncontrol <- as.numeric(phen$V4[24])-ncase
    K <- P <- ncase/(ncase+ncontrol)
    res$lder$h2.lia <- cvt.h2(K,P,res$lder$h2)
    res$lder$h2.sd.lia <- cvt.h2(K,P,res$lder$h2.sd)
    res$ldsc$h2.lia <- cvt.h2(K,P,res$ldsc$h2)
    res$ldsc$h2.sd.lia <- cvt.h2(K,P,res$ldsc$h2.sd)
    print(res)
    save(res,file='res_new_lder.lia.RData')
  }
  for(trait in c('T2D','ASTHMA','CAD','SCZ')){
    setwd('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/')
    dir.create(trait)
    print(trait)
    setwd(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait))
    load('res_new_lder.lia.RData')
    phen <- fread('plink.log',fill=T)
    ncase <- as.numeric(phen$V4[25])
    ncontrol <- as.numeric(phen$V4[24])-ncase
    K <- P <- ncase/(ncase+ncontrol)
    res$lder$h2 <- res$h2
    res$lder$h2.sd <- res$h2.sd
    res$lder$h2.lia <- cvt.h2(K,P,res$lder$h2)
    res$lder$h2.sd.lia <- cvt.h2(K,P,res$lder$h2.sd)
    print(res)
    save(res,file='res_new_lder.lia.RData')
  }
  
}


get.res <- function(trait,binary=F,lia=T){
  if(!binary){
    load(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait,'/res_new_lder.RData'))
    return(list(lder=res$h2,
                lder.sd=res$h2.sd
                ))
  }else{
    load(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait,'/res_new_lder.lia.RData'))
    if(lia){return(list(lder=res$lder$h2.lia,
                        lder.sd=res$lder$h2.sd.lia))
    }else{  return(list(lder=res$lder$h2,
                        lder.sd=res$lder$h2.sd,))     
    }
  }
}
k <- 1
for(trait in c('LDL','HDL','TG','TC','HEIGHT','BMI')){
  if(k==1){
    dat <- unlist(get.res(trait))
    k <- k+1
  }else{
    dat <- rbind(dat,unlist(get.res(trait)))
  }
}
setwd('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits')
write.csv(dat,'continuous_newinf.csv',quote=F,row.names=F)

k <- 1
for(trait in c('ASTHMA','CAD','SCZ','T2D')){
  if(k==1){
    dat <- unlist(get.res(trait,binary=T))
    k <- k+1
  }else{
    dat <- rbind(dat,unlist(get.res(trait,binary=T)))
  }
}
setwd('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits')
write.csv(dat,'binary_newinf.csv',quote=F,row.names=F)

