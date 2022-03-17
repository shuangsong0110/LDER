source('/home/songs/genomic_control/ukb/20210204/ldsc_new.R')
path0 <- '/home/songs/genomic_control/ukb/analysis/ukb_ld/'
setwd(path0)
library(data.table)
bim <- fread('/home/songs/genomic_control/ukb/dat/alldat.qc.bim')
fam <- fread('/home/songs/genomic_control/ukb/dat/alldat.qc.fam')

source('/home/songs/genomic_control/20210714_code.R')

write.assoc <- function(n){
  set.seed(0)
  path <- paste0(path0,'n',n,'/')
  dir.create(path)
  setwd(path)
  write.table(fam[sample(nrow(fam),n),1:2],'sample.txt',quote=F,row.names=F,col.names=F)
  if(!file.exists('dat00.bed')){
  system(paste0('/home/songs/plink1.9/plink --bfile  /home/songs/genomic_control/ukb/dat/alldat.qc ',
                '--keep ',path,'sample.txt --make-bed --out ',path,'dat00'))}
  for(times in 1:50){

    for(alpha in c(0.01,0.005)){
      path <- paste0(path0,'n',n,'/')
      dir.create(path)
      setwd(path)
   
      bim00 <- fread('dat00.bim')
      m <- nrow(bim00)
      set.seed(times)
      beta <- rep(0,m)
      ind <- sample(m,alpha*m)
      write.table(bim00$V2[ind],'snp.txt',quote=F,row.names=F,col.names=F)
      system(paste0('/home/songs/plink1.9/plink --bfile  ',path,'dat00 ',
                    '--extract ',path,'snp.txt --make-bed --out ',path,'dat001'))
      library(EBPRS)
      dat <- read_plink(paste0(path,'dat001'))
      set.seed(times)
      X <- dat$bed
      X[which(is.na(X))] <- 0
      X <- scale(X)
      X[which(is.na(X))] <- 0
      h2 <- 0.5
      m <- nrow(dat$bim)
      beta <- rnorm(length(ind),0,sqrt(h2/m))
      e <- rnorm(n,0,sqrt(1-h2))
      y <- X%*%beta+e
      dat$fam$V6 <- scale(y)
      write.table(dat$fam,'dat00.fam',quote=F,row.names=F,col.names=F)
      dir.create(paste0('./h2_',h2,'/alpha_',alpha),recursive = T)
      system(paste0('/home/songs/plink1.9/plink --bfile ',path,'dat00 --assoc --out ',path,'h2_',h2,'/alpha_',alpha,
                    '/times',times))
      ####
     }
 # }
  }
}
#write.assoc(2000,1)
library(parallel)
mclapply(c(2000,5000,1e4,2e4,5e4,1e5),write.assoc,mc.cores=6)

h2 <- 0.5
run.result <- function(n,alpha,times,ld.type='ukb',twostage=T){
  set.seed(0)
  path <- paste0(path0,'n',n,'/')
  dir.create(path)
  setwd(path)
  assoc <- fread(paste0(path,'h2_',h2,'/alpha_',alpha,
                 '/times',times,'.qassoc'))
  
  bim00 <- fread('dat00.bim')
  bim01 <- bim00[,c(2,5,6)]
  colnames(bim01) <- c('SNP','a0','a1')
  assoc <- merge(assoc,bim01)
  assoc$z <- -sign(assoc$BETA)*qnorm(assoc$P/2)
  colnames(assoc)[which(colnames(assoc)=='SNP')] <- 'snp'
  colnames(assoc)[which(colnames(assoc)=='CHR')] <- 'chr'
  colnames(assoc)[which(colnames(assoc)=='BP')] <- 'bp'
  colnames(assoc)[which(colnames(assoc)=='NMISS')] <- 'N'
  if(ld.type=='1kg'){
    print('onekg')
  ldpath <- '/home/songs/genomic_control/ukb/20210204/1kg_sub/LD/'
  ldpath.shrink <-  '/home/songs/genomic_control/ukb/20210204/1kg_sub/LD.shrink/'
  res <- ld.all(assoc,assoc$N[1],ldpath=ldpath,ldpath.shrink=ldpath.shrink,method='both',
                type='boot',rough=T,twostage=twostage,cores=10,n.ld=489,n.bs=500)
  }else if(ld.type=='ukb'){
    ldpath <- '/home/songs/genomic_control/ukb/20210204/ukb/LD/'
    ldpath.shrink <-  '/home/songs/genomic_control/ukb/20210204/ukb/LD.shrink/'
    res <- ld.all(assoc,assoc$N[1],ldpath=ldpath,ldpath.shrink=ldpath.shrink,method='both',
                  type='boot',rough=F,twostage=twostage,cores=10,n.bs=500)
  }
  
  return(res)
  
}
if(F){
system.time(
tt <- run.result(2000,0.005,3,ld.type='ukb',twostage=T)
)
}

library(parallel)
h2.all <- method <- nn <- sd.all <- c()
alpha <- 0.005
### ukb
for (n in c(2000,5000,1e4,2e4,5e4)){
  res.ukb <- mclapply(1:50,run.result,n=n,alpha=alpha,ld.type='ukb',mc.cores=20)
  lder11 <- ldsc11 <- lder11.sd <- ldsc11.sd <- c()
  for(k in 1:50){
    lder11[k] <- res.ukb[[k]]$lder$h2
    ldsc11[k] <- res.ukb[[k]]$ldsc$h2
    lder11.sd[k] <- res.ukb[[k]]$lder$h2.sd
    ldsc11.sd[k] <- res.ukb[[k]]$ldsc$h2.sd
  }
  nn <- c(nn,rep(n,100))
  h2.all <- c(h2.all,lder11,ldsc11)
  sd.all <- c(sd.all,lder11.sd,ldsc11.sd)
  method <- c(method,rep('LDER',50),rep('LDSC',50))
  dat.ukb <- data.frame(Method=method,h2=h2.all,n=nn,sd=sd.all)
  setwd('/home/songs/genomic_control/ukb/analysis/ukb_ld/plot_data_real_simu/')
  save(nn,h2.all,sd.all,method,dat.ukb,file=paste0('dat.ukb.h2_',h2,'_alpha_',alpha,
                                     '.1120.RData'))
}

### 1kg
h2.all <- method <- nn <- sd.all <- c()
for (n in c(2000,5000,1e4,2e4,5e4)){
  res.ukb <- mclapply(1:50,run.result,n=n,alpha=alpha,ld.type='1kg',mc.cores=10)
  lder11 <- ldsc11 <- lder11.sd <- ldsc11.sd <- c()
  for(k in 1:50){
    lder11[k] <- res.ukb[[k]]$lder$h2
    ldsc11[k] <- res.ukb[[k]]$ldsc$h2
    lder11.sd[k] <- res.ukb[[k]]$lder$h2.sd
    ldsc11.sd[k] <- res.ukb[[k]]$ldsc$h2.sd
  }
  nn <- c(nn,rep(n,100))
  h2.all <- c(h2.all,lder11,ldsc11)
  sd.all <- c(sd.all,lder11.sd,ldsc11.sd)
  method <- c(method,rep('LDER',50),rep('LDSC',50))
  dat.onekg <- data.frame(Method=method,h2=h2.all,n=nn,sd=sd.all)
  setwd('/home/songs/genomic_control/ukb/analysis/ukb_ld/plot_data_real_simu/')
  save(nn,h2.all,sd.all,dat.onekg,file=paste0('dat.onekg.h2_',h2,'_alpha_',alpha,
                           '.1120.RData'))
}

h2 <- 0.5
alpha <- 0.005
dir.create('/home/songs/genomic_control/ukb/analysis/ukb_ld/plot_data_real_simu/')
setwd('/home/songs/genomic_control/ukb/analysis/ukb_ld/plot_data_real_simu/')
load(paste0('dat.onekg.h2_',h2,'_alpha_',alpha,
            '.1120.RData'))
load(paste0('dat.ukb.h2_',h2,'_alpha_',alpha,
            '.1120.RData'))
save(dat.ukb,dat.onekg,file=paste0('h2_',h2,'_alpha_',alpha,
                                   '.1120.RData'))


#### summary
h2 <- 0.5
alpha=0.005
setwd('/home/songs/genomic_control/ukb/analysis/ukb_ld/plot_data_real_simu/')
load(paste0('h2_',h2,'_alpha_',alpha,
            '.1120.RData'))
### precision
dat <- dat.ukb
dat$prec <- 1/dat$sd
for(n in c(2000,5000,1e4,2e4,5e4)){
  #tt <-sqrt( (dat$h2[which((dat$Method=='LDSC')&(dat$n==n))]-0.5)^2)
  tt <- dat$prec[which((dat$Method=='LDER')&(dat$n==n))]
  print(n)
  print(mean(tt))
  print(sd(tt))
}

### used sd instead of se

dat <- dat.onekg
dat$prec <- 1/dat$sd
for(n in c(5000,1e4,2e4,5e4)){
  #tt <-sqrt( (dat$h2[which((dat$Method=='LDSC')&(dat$n==n))]-0.5)^2)
  tt <- 1/sd(dat$h2[which((dat$Method=='LDER')&(dat$n==n))])
  print(n)
  print(tt)
}


### rmse

h2 <- 0.5
alpha=0.005
setwd('/home/songs/genomic_control/ukb/analysis/ukb_ld/plot_data_real_simu/')
load(paste0('h2_',h2,'_alpha_',alpha,
            '.1120.RData'))
dat <- dat.ukb
for(n in c(2000,5000,1e4,2e4,5e4)){
  tt <-( (dat$h2[which((dat$Method=='LDER')&(dat$n==n))]-0.5)^2)
  print(n)
  print(sqrt(mean(tt)))
}

for(n in c(2000,5000,1e4,2e4,5e4)){
  tt <-( (dat$h2[which((dat$Method=='LDSC')&(dat$n==n))]-0.5)^2)
  print(n)
  print(sqrt(mean(tt)))
}



############################

get.p <- function(dat,n){
  tt <- dat[dat$n==n,]
  ts <- var.test(h2~Method,tt,alternative='less')
  return(ts$p.value)
}
get.pv <- Vectorize(get.p,'n')


get.plot.pv <- function(h2,alpha){
  library(ggpubr)
  library(rstatix)
  library(ggsci)
  load(paste0('./plot_data_real_simu/h2_',h2,'_alpha_',alpha,
              '.plot_dat50.new_2stage.RData'))
  dat <- dat.ukb
  library(ggpubr)
  library(rstatix)
  dat$N1 <- paste0('n = ',dat$n)
  dat$N1 <- factor(dat$N1,levels=unique(dat$N1))
  
  bxp <- ggboxplot(
    dat, x = "Method", y = "h2", fill = "Method", 
    facet.by = "N1",scales="free"
  )
  stat.test <- dat %>%
    group_by(N1) %>%
    t_test(h2 ~ Method) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()
  stat.test  ### this test is not right... should change to variance test
  pvalue <- get.pv(dat,unique(dat$n))
  pvalue1 <- format(pvalue,scientific = T,digits = 1)
  stat.test$p.adj <- pvalue1
  
  library(ggsci)
  # Make facet and add p-values
  stat.test <- stat.test %>% add_xy_position(x = "Method")
  p2 <- bxp +   
    stat_pvalue_manual(
      stat.test, bracket.nudge.y = -1, hide.ns = TRUE,
      label = "p={p.adj}"
    )+
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
    scale_fill_npg()+
    xlab(' ')+
    ylab(expression(Estimated~h^2))+
  theme(panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_line(size=0.5),
        legend.title =  element_text(size=13),
        legend.text =  element_text(size=12),
        axis.title=element_text(size=13,face="bold"),
        legend.position = 'bottom')+
    geom_segment(x = 0, xend = 4, y = h2, yend = h2,col='#a58acc',size=.5,lty=2)
  return(p2)
}

pp <- get.plot.pv(0.5,0.01)
pp
pdf(paste0('./figures/real_geno_ukb_h2_',h2,'_alpha_',alpha,'.pdf'),height=5.5,width = 8,onefile=F)
get.plot.pv(0.5,0.01)
dev.off()

library(ggsci)
#'{}^0*italic(D)~plain((rarefied))'
get.plot.both <- function(h2,alpha){
  load(paste0('./plot_data_real_simu/h2_',h2,'_alpha_',alpha,
              '.plot_dat50.new_2stage.RData'))
  dat <- dat.ukb
  if(T){
    dat <- dat[dat$n!=2000,]
  }
  dat$N1 <- paste0('n = ',dat$n)
  dat$N1 <- factor(dat$N1,levels=unique(dat$N1))
  p1 <-  ggplot(dat, aes(x = Method, y = h2)) + 
    geom_boxplot(aes(fill = Method),
                 position = position_dodge(0.9)) +
    scale_fill_npg()+scale_color_npg()+theme_bw()+
    #theme(legend.position = 'none')+ 
    xlab(' ')+
    ylab(expression(Estimated~h^2))+theme_bw()+
    theme(panel.grid=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line=element_line(size=0.5),
          legend.title =  element_text(size=13),
          legend.text =  element_text(size=12),
          axis.title=element_text(size=13,face="bold"))+
    # facet_wrap(~N, labeller=label_parsed,nrow=1)+
    geom_segment(x = 0, xend = 4, y = h2, yend = h2,col='#a58acc',size=.5,lty=2)+
    #facet_wrap(~N1,scales="free",nrow=1) +
    facet_wrap(~N1,nrow=1) +
    stat_summary(fun=mean,
                 colour="darkred", geom="point", 
                 shape=18, size=3,show.legend = FALSE)
  if(pv){
    pvalue <- get.pv(dat,unique(dat$n))
    yy <- dat$h2[dat$Method==dat$Method[1]]+0.1
    xx1 <- 1-0.15
    xx2 <- 1+0.15
  library(ggsignif)
  p11 <- p1+geom_signif(y_position = yy, xmin =xx1, xmax = xx2,
                         annotations = paste0('p = ',pvalue), tip_length = 0.01)
  }
  p11
  dat <- dat.onekg
  dat$N1 <- paste0('n = ',dat$n)
  dat$N1 <- factor(dat$N1,levels=unique(dat$N1))
  p2 <-  ggplot(dat, aes(x = Method, y = h2)) + 
    geom_boxplot(aes(fill = Method),
                 position = position_dodge(0.9)) +
    scale_fill_npg()+scale_color_npg()+theme_bw()+
    #theme(legend.position = 'none')+ 
    xlab(' ')+
    ylab(expression(Estimated~h^2~LD[1000][G]))+theme_bw()+
    theme(panel.grid=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line=element_line(size=0.5),
          legend.title =  element_text(size=13),
          legend.text =  element_text(size=12),
          axis.title=element_text(size=13,face="bold"))+
    # facet_wrap(~N, labeller=label_parsed,nrow=1)+
    geom_segment(x = 0, xend = 4, y = h2, yend = h2,col='#a58acc',size=.5,lty=2)+
    #facet_wrap(~N1,scales="free",nrow=1) +
    facet_wrap(~N1,nrow=1) +
    stat_summary(fun=mean,
                 colour="darkred", geom="point", 
                 shape=18, size=3,show.legend = FALSE)
  p2
  p3 <- lemon::grid_arrange_shared_legend(p1,p2,nrow=2,ncol=1)
  p3
}
library(ggplot2)
p3 <- get.plot.both(h2=0.5,alpha=0.01)
p3

h2 <- 0.5;alpha <- 0.01
pdf(paste0('./figures/real_geno_h2_',h2,'_alpha_',alpha,'.1122.pdf'),height=6,width = 14,onefile=F)
get.plot.both(h2=h2,alpha=alpha)
dev.off()







