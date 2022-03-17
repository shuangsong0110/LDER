source('/home/songs/genomic_control/ukb/20210204/ldsc_new.R')

path0 <- '/home/songs/genomic_control/20201218simu/new/revision220127_fixR/'
dir.create(path0)
setwd(path0)
### save statistics
set.seed(0)
m <- 1e5
m.bl <- m/100
m0 <- m/m.bl
p.all <- runif(m.bl,0.1,0.9)
## get result for 

source('/home/songs/genomic_control/ukb/20210204/ldsc_new.R')

path0 <- '/home/songs/genomic_control/20201218simu/new/revision220127_fixR/'
setwd(path0)


get.res <- function(N,m,alpha=0.05,h2=0.5,inf=1,rough=F,twostage=T,a=NULL){
  a0 <- a
  path <- paste0(path0,'/inf_',inf,'/h2_',h2,'/alpha_',alpha,'/m_',m,'/N_',N)
  n.gwas=N
  res.all <- matrix(0,50,4)
  load(paste0(path,'/stats.RData'))
  for(times in 1:50){
    res <- c()
    x1 <- x.all[[times]] 
    z1 <- z.all[[times]]
    lam1 <- lam.all[[times]] 
    ldsc1 <- ld.all[[times]] 
    #m <- length(x1)
    newGC <-  1+n.gwas*calH2.new1(x1, lam1,   N=n.gwas,a=a0, rough=F)$a
    #print(newGC)
    if(newGC<1.05){
      s2thld <- 3.5*sqrt(m/n.gwas)+2.5 
    }else{
      s2thld <- quantile(x1^2, 1-5/m)
    }
    idx <- (x1^2<=s2thld)
    temp1 <-  calH2.new1(x1[idx], lam1[idx],   N=n.gwas,a=a0, rough=F)$a
    a <- max(temp1,0)
    a <- max(calH2.new1(x1, lam1,   N=n.gwas,a=a0, rough=F)$a,0)
    if(twostage){
      h2 <- calH2.new1(x1, lam1,   N=n.gwas, a=max(temp1,0), rough=rough)$h2
    }else{
      h2 <- calH2.new1(x1, lam1,   N=n.gwas, a=a0, rough=rough)$h2
      a <- calH2.new1(x1, lam1,   N=n.gwas, a=a0, rough=rough)$a
    }
    
    res[1] <- h2
    ## 2. two-stage
    #temp <-  calH2.new1(x1[(x1^2)<s2thld], lam1[(x1^2)<s2thld],   N=n,a=NULL)$a
    res[2] <-a*n.gwas+1
    s1thld <- 30
    if(twostage){
      a <-  calH2(z1[(z1^2)<s1thld], ldsc1[(z1^2)<s1thld],   N=n.gwas,a=a0)$a
      h2 <- calH2(z1, ldsc1,   N=n.gwas,a=a)$h2
    }else{
      a <- calH2(z1, ldsc1,   N=n.gwas,a=a0)$a
      h2 <-  calH2(z1, ldsc1,   N=n.gwas,a=a0)$h2
    }
    
    res[3] <- h2
    ## 1. origin
    res[4] <-a*n.gwas+1
    res.all[times,] <- res
  }
  res.all <- data.frame(res.all)
  colnames(res.all) <- c('lder.h2','lder.a','ldsc.h2',
                         'ldsc.a')
  return(res.all)
}

res1 <- get.res(5e3,1e5,alpha=0.01,inf=1,h2=0.05,twostage=F)  
res2 <- get.res(5e4,1e5,alpha=0.01,inf=1,h2=0.5,twostage=F)  

#apply(res,2,mean)
apply(res,2,sd)


path0 <- '/home/songs/genomic_control/20201218simu/new/revision220127_fixR/'
setwd(path0)
for(h2 in c(0.2,0.5,0.05,0.1)){
  for(alpha in c(0.005,0.01,0.05,0.1)){
    for(m in c(1e5)){
      #  for(N in c(2000,5000,10000,2e4,5e4,1e5,2e5)){
      for (inf in c(1,1.1)){
        Method <- rep(rep(c('LDER','LDSC'),each=50),4)
        N <- rep(c(5000,1e4,2e4,5e4),each=100)
        h2.est <- c()
        a.est <- c()
        print(paste0(path0,'/plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
                     '_inf_',inf,'.plot_dat50.new_1stage_inf1stage.RData'))
        for(n.gwas in c(5000,1e4,2e4,5e4)){
          res.all <- get.res(n.gwas,m,alpha=alpha,inf=inf,h2=h2,twostage = F)  
          h2.est <- c(h2.est,res.all[,1],res.all[,3])
          #res.all <- get.a(n.gwas,m,alpha=alpha,inf=inf)  
          a.est <- c(a.est,res.all[,2],res.all[,4])
        }
        dat <- data.frame(Method,N,h2.est,a.est)
        dir.create(paste0(path0,'/plot_dat'))
        save(dat,file=paste0(path0,'/plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
                             '_inf_',inf,'.plot_dat50.new_1stage_inf1stage.RData'))
        #get.stats(N=N,m=m,alpha=alpha,h2=h2,times.all=20,inf=inf)
        #}
      }
    }
  }
}

#### plot
library(ggplot2)
library(ggsci)
#'{}^0*italic(D)~plain((rarefied))'
if(F){
  get.plot <- function(h2,alpha,m,inf,type){
    load(paste0('./0208plot_data/h2_',h2,'_alpha_',alpha,'_m_',m,
                '_inf_',inf,'.plot_dat50.new_1stage_inf1stage.RData'))
    dat$N1 <- paste0('n = ',dat$N)
    dat$N1 <- factor(dat$N1,levels=unique(dat$N1))
    if(type=='h2'){
      p <- ggplot(dat, aes(x = Method, y = h2.est)) + 
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
              legend.title =  element_text(size=14),
              legend.text =  element_text(size=12),
              axis.title=element_text(size=15,face="bold"))+
        # facet_wrap(~N, labeller=label_parsed,nrow=1)+
        geom_segment(x = 0, xend = 4, y = h2, yend = h2,col='#a58acc',size=.5,lty=2)+
        facet_wrap(~N1,scales="free") +
        stat_summary(fun=mean,
                     colour="darkred", geom="point", 
                     shape=18, size=3,show.legend = FALSE)
      
    }
    if(type=='inf'){
      p <- ggplot(dat, aes(x = Method, y = a.est)) + 
        geom_boxplot(aes(fill = Method),
                     position = position_dodge(0.9)) +
        scale_fill_npg()+scale_color_npg()+theme_bw()+
        #theme(legend.position = 'none')+ 
        xlab(' ')+
        ylab(expression(Estimated~inflation))+theme_bw()+
        theme(panel.grid=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.line=element_line(size=0.5),
              legend.title =  element_text(size=14),
              legend.text =  element_text(size=12),
              axis.title=element_text(size=15,face="bold"))+
        # facet_wrap(~N, labeller=label_parsed,nrow=1)+
        geom_segment(x = 0, xend = 4, y = inf, yend = inf,col='#a58acc',size=.5,lty=2)+
        facet_wrap(~N1,scales="free") +
        stat_summary(fun=mean,
                     colour="darkred", geom="point", 
                     shape=18, size=3,show.legend = FALSE)
      
    }
    return(p)
  }
}
get.plot.both <- function(h2,alpha,m,inf){
  library(scales)
  load(paste0('./plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
              '_inf_',inf,'.plot_dat50.new_1stage_inf1stage.RData'))
  # dat <- dat[dat$N!=1e5,] ## not include 1e5 here
  # dat <- dat[dat$N!=2000,] ## not include 1e5 here
  dat$N1 <- paste0('n = ',dat$N)
  dat$N1 <- factor(dat$N1,levels=unique(dat$N1))
  p1 <-  ggplot(dat, aes(x = Method, y = h2.est)) + 
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
          axis.title=element_text(size=13,face="bold"),
          # The new stuff
          strip.text = element_text(size = 14))+
    # facet_wrap(~N, labeller=label_parsed,nrow=1)+
    geom_segment(x = 0, xend = 4, y = h2, yend = h2,col='#a58acc',size=.5,lty=2)+
    #facet_wrap(~N1,scales="free",nrow=1) +
    facet_wrap(~N1,nrow=1) +
    stat_summary(fun=mean,
                 colour="darkred", geom="point", 
                 shape=18, size=3,show.legend = FALSE)+
    scale_y_continuous(breaks= pretty_breaks())
  p1
  p2 <- ggplot(dat, aes(x = Method, y = a.est)) + 
    geom_boxplot(aes(fill = Method),
                 position = position_dodge(0.9)) +
    scale_fill_npg()+scale_color_npg()+theme_bw()+
    #theme(legend.position = 'none')+ 
    xlab(' ')+
    ylab(expression(Estimated~inflation))+theme_bw()+
    theme(panel.grid=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line=element_blank(),
          #axis.line=element_line(size=0.5),
          legend.title =  element_text(size=13),
          legend.text =  element_text(size=12),
          axis.title=element_text(size=13,face="bold"),
          # The new stuff
          strip.text = element_text(size = 14))+
    # facet_wrap(~N, labeller=label_parsed,nrow=1)+
    geom_segment(x = 0, xend = 4, y = inf, yend = inf,col='#a58acc',size=.5,lty=2)+
    #facet_wrap(~N1,scales="free",nrow=1) +
    facet_wrap(~N1,nrow=1) +
    stat_summary(fun=mean,
                 colour="darkred", geom="point", 
                 shape=18, size=3,show.legend = FALSE)+
    scale_y_continuous(breaks= pretty_breaks())
  p3 <- lemon::grid_arrange_shared_legend(p1,p2,nrow=2,ncol=1)
  return(p3)
}




# p1 <- get.plot(h2=0.5,alpha=0.01,m=1e5,inf=1,type='h2')
# p1
# p2 <- get.plot(h2=0.5,alpha=0.01,m=1e5,inf=1,type='inf')
# p2

p3 <- get.plot.both(h2=0.5,alpha=0.05,m=1e5,inf=1)
p3
library(ggplot2)
for(h2 in c(0.2,0.5)){
  for(alpha in c(0.005,0.01)){
    m <- 1e5
    for(inf in c(1,1.1)){
      print(paste0('./figures/h2_',h2,'_alpha_',alpha,'_m_',m,
                   '_inf_',inf,'_inf1st.pdf'))
      pdf(paste0('./figures/h2_',h2,'_alpha_',alpha,'_m_',m,
                 '_inf_',inf,'_inf1st.new.pdf'),height=5.5,width = 9,onefile = F)
      get.plot.both(h2=h2,alpha=alpha,m=m,inf=inf)
      dev.off()
    }}}



get.plot.both.compare <- function(h2,alpha,m,inf){
  library(scales)
  load(paste0('./plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
              '_inf_',inf,'.plot_dat50.new_1stage_inf1stage.RData'))
  head(dat)
  # dat <- dat[dat$N!=1e5,] ## not include 1e5 here
  # dat <- dat[dat$N!=2000,] ## not include 1e5 here
  dat$N1 <- paste0('n = ',dat$N)
  dat$N1 <- factor(dat$N1,levels=unique(dat$N1))
  dat.1st <- dat
  dat.1st$Method <- paste0(dat.1st$Method ,'(1-stage)')
  
  load(paste0('./plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
              '_inf_',inf,'.plot_dat50.new_2stage_inf1stage.all.RData'))
  # dat <- dat[dat$N!=1e5,] ## not include 1e5 here
  # dat <- dat[dat$N!=2000,] ## not include 1e5 here
  dat$N1 <- paste0('n = ',dat$N)
  dat$N1 <- factor(dat$N1,levels=unique(dat$N1))
  dat <- dat[dat$Method%in%c('LDER','LDSC'),]
  dat$Method <- paste0(dat$Method ,'(2-stage)')
  dat <- rbind(dat,dat.1st)
  dat$Method <- factor(dat$Method,levels=unique(dat$Method)[c(1,3,2,4)])
  
  p1 <-  ggplot(dat, aes(x = Method, y = h2.est)) + 
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
          axis.title=element_text(size=13,face="bold"),
          # The new stuff
          strip.text = element_text(size = 14))+
    # facet_wrap(~N, labeller=label_parsed,nrow=1)+
    geom_segment(x = 0, xend = 4, y = h2, yend = h2,col='#a58acc',size=.5,lty=2)+
    #facet_wrap(~N1,scales="free",nrow=1) +
    facet_wrap(~N1,nrow=1) +
    stat_summary(fun=mean,
                 colour="darkred", geom="point", 
                 shape=18, size=3,show.legend = FALSE)+
    scale_y_continuous(breaks= pretty_breaks())
  p1
  p2 <- ggplot(dat, aes(x = Method, y = a.est)) + 
    geom_boxplot(aes(fill = Method),
                 position = position_dodge(0.9)) +
    scale_fill_npg()+scale_color_npg()+theme_bw()+
    #theme(legend.position = 'none')+ 
    xlab(' ')+
    ylab(expression(Estimated~inflation))+theme_bw()+
    theme(panel.grid=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line=element_blank(),
          #axis.line=element_line(size=0.5),
          legend.title =  element_text(size=13),
          legend.text =  element_text(size=12),
          axis.title=element_text(size=13,face="bold"),
          # The new stuff
          strip.text = element_text(size = 14))+
    # facet_wrap(~N, labeller=label_parsed,nrow=1)+
    geom_segment(x = 0, xend = 4, y = inf, yend = inf,col='#a58acc',size=.5,lty=2)+
    #facet_wrap(~N1,scales="free",nrow=1) +
    facet_wrap(~N1,nrow=1) +
    stat_summary(fun=mean,
                 colour="darkred", geom="point", 
                 shape=18, size=3,show.legend = FALSE)+
    scale_y_continuous(breaks= pretty_breaks())
  p3 <- lemon::grid_arrange_shared_legend(p1,p2,nrow=2,ncol=1)
  return(p3)
}
p3 <- get.plot.both.compare(h2=0.5,alpha=0.05,m=1e5,inf=1.1)
p3
pdf('simu_fixR_1stage1.pdf',height=5.5,width = 9,onefile = F)
get.plot.both.compare(h2=0.5,alpha=0.05,m=1e5,inf=1)
dev.off()
pdf('simu_fixR_1stage2.pdf',height=5.5,width = 9,onefile = F)
get.plot.both.compare(h2=0.5,alpha=0.05,m=1e5,inf=1.1)
dev.off()




