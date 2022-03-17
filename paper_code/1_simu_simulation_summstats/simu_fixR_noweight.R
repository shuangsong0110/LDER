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
save(p.all,file=paste0('p.all.m',m,'.RData'))

get.stats.times<- function(N,m,alpha=0.05,h2=0.5,times,inf=1){
  print(m)
  path <- paste0(path0,'/inf_',inf,'/h2_',h2,'/alpha_',alpha,'/m_',m,'/N_',N)
  dir.create(path,recursive=T)
  load(paste0(path0,'/p.all.m',m,'.RData'))
  if(!file.exists(paste0(path,'/stats.RData'))){
    print(N)
    res1 <- res2 <- res1.0 <- res2.0 <- c()
    x.all <- z.all <- lam.all <- ld.all <- list()
    #  for(alpha in  c(0.05)){
    #  for(times in times.all){
    set.seed(times)
    m.bl <- m/100
    library(Matrix)
    #temp <- 1000
    beta <- rep(0,m)
    ind <- sample(m,alpha*m)
    beta[ind] <- rnorm(length(ind),0,sqrt(h2/length(ind)))
    z <- ev <- x <- ldsc <- c()
    for(i in 1:m.bl){
      #print(i)
      # m0 <- m/m.bl
      # p <- runif(1,0.1,0.9)
      p <- p.all[i]
      R <- p^(abs(outer(1:m0, 1:m0, "-")))
      library(mvtnorm)
      beta0 <- beta[((i-1)*m0+1):(i*m0)]
      z.mean <- sqrt(N) * R %*% beta0
      z0 <- rmvnorm(1, mean=z.mean, sigma=inf * R)
      temp <- eigen(R)
      U <- temp$vectors
      V <- diag(temp$values)
      V.inv <- diag(1/temp$values)
      x0 <- sqrt(V.inv)%*%t(U)%*%t(z0)
      ldsc0 <- apply(R^2, 2, sum)
      
      z <- c(z,z0)
      x <- c(x,x0)
      ev <- c(ev,temp$values)
      ldsc <- c(ldsc,ldsc0)
    }
    s1thld <- 30
    ### eigen.shrink
    x1 <- x;lam1 <- ev;z1 <- z;ld1 <- ldsc;
    return(list(x1=x1,z1=z1,lam1=lam1,ld1=ld1))
  }
}
get.stats.all <-  function(N,m,alpha=0.05,h2=0.5,times.all=1:20,inf=1,cores=30){
  print(m)
  path <- paste0(path0,'/inf_',inf,'/h2_',h2,'/alpha_',alpha,'/m_',m,'/N_',N)
  dir.create(path,recursive=T)
  
  library(parallel)
  res <- mclapply(times.all,get.stats.times,N=N,m=m,alpha=alpha,h2=h2,inf=inf,mc.cores=cores)
  res1 <- res2 <- res1.0 <- res2.0 <- c()
  x.all <- z.all <- lam.all <- ld.all <- list()
  
  for(j in 1:length(times.all)){
    x.all[[j]] <- res[[j]]$x1
    z.all[[j]] <- res[[j]]$z1
    lam.all[[j]] <- res[[j]]$lam1
    ld.all[[j]] <- res[[j]]$ld1
  }
  save(x.all,z.all,lam.all,ld.all,file=paste0(path,'/stats.RData'))
  
}
get.stats.all(1e5,1e5,alpha,times.all=1)
## new
for(h2 in c(0.05,0.1,0.2,0.5)){
  for(alpha in c(0.005,0.01,0.05,0.1)){
    for(m in c(1e5)){
      for (inf in c(1,1.1)){
        for(N in c(5000,10000,2e4,5e4)){
          get.stats.all(N,m,alpha=alpha,h2=h2,times.all=1:50,inf=inf,cores=26)
        }
      }
    }
  }
}

## get result for 

source('/home/songs/genomic_control/ukb/20210204/ldsc_new.R')

path0 <- '/home/songs/genomic_control/20201218simu/new/revision220127_fixR/'
setwd(path0)


calH2.new1.noweight <-function(z, ldsc, N, a=NULL, rough=F){ ## ldsc
  if(length(z)!=length(ldsc)){
    stop('ldsr.calH2:lengths of z and ldsc are not matched!')
  }
  M <- length(z)
  chi2 <- z^2
  updateFunc <- function(coeff, w){
    temp1 <- min(max(0,coeff[1]),N/M)
    #temp1 <- coeff[1]
    w1 <- 1/(1+temp1*ldsc)^2
    if(rough){
      w2 <- 1/(1+exp(-5*(ldsc-1)))
    }else{
      w2 <- pmin(1, ldsc)
    }
    return(rep(1,M))
  }
  if(!is.null(a)){
    int <- N*a+1.
    y <- chi2-int
    coeff <- irwls(as.matrix(ldsc), y, updateFunc)
  }else{
    x <- cbind(ldsc, rep(1, length(ldsc)))
    y <- chi2
    coeff <- irwls(x, y, updateFunc)
    a <- (coeff[[2]]-1)/N
  }
  h2est <- coeff[1]*M/N
  #h2est <- min(max(0,coeff[1]),N/M)*M/N
  return(list(h2=h2est, a=a))
} 

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
    newGC <-  1+n.gwas*calH2.new1.noweight(x1, lam1,   N=n.gwas,a=a0, rough=F)$a
    #print(newGC)
    if(newGC<1.05){
      s2thld <- 3.5*sqrt(m/n.gwas)+2.5 
    }else{
      s2thld <- quantile(x1^2, 1-5/m)
    }
    idx <- (x1^2<=s2thld)
    temp1 <-  calH2.new1.noweight(x1[idx], lam1[idx],   N=n.gwas,a=a0, rough=F)$a
    a <- max(temp1,0)
    a <- max(calH2.new1.noweight(x1, lam1,   N=n.gwas,a=a0, rough=F)$a,0)
    if(twostage){
      h2 <- calH2.new1.noweight(x1, lam1,   N=n.gwas, a=max(temp1,0), rough=rough)$h2
    }else{
      h2 <- calH2.new1.noweight(x1, lam1,   N=n.gwas, a=a0, rough=rough)$h2
      a <- calH2.new1.noweight(x1, lam1,   N=n.gwas, a=a0, rough=rough)$a
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

res1 <- get.res(5e3,1e5,alpha=0.01,inf=1,h2=0.05,twostage=T)  
res2 <- get.res(5e4,1e5,alpha=0.01,inf=1,h2=0.5,twostage=T)  



path0 <- '/home/songs/genomic_control/20201218simu/new/revision220127_fixR/'
setwd(path0)
for(h2 in c(0.5)){
  for(alpha in c(0.05)){
    for(m in c(1e5)){
      #  for(N in c(2000,5000,10000,2e4,5e4,1e5,2e5)){
      for (inf in c(1,1.1)){
        Method <- rep(rep(c('LDER','LDSC'),each=50),4)
        N <- rep(c(5000,1e4,2e4,5e4),each=100)
        h2.est <- c()
        a.est <- c()
        print(paste0(path0,'/plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
                     '_inf_',inf,'.plot_dat50.new_2stage_inf1stage.RData'))
        for(n.gwas in c(5000,1e4,2e4,5e4)){
          res.all <- get.res(n.gwas,m,alpha=alpha,inf=inf,h2=h2)  
          h2.est <- c(h2.est,res.all[,1],res.all[,3])
          #res.all <- get.a(n.gwas,m,alpha=alpha,inf=inf)  
          a.est <- c(a.est,res.all[,2],res.all[,4])
        }
        dat <- data.frame(Method,N,h2.est,a.est)
        dir.create(paste0(path0,'/plot_dat'))
        save(dat,file=paste0(path0,'/plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
                             '_inf_',inf,'.plot_dat50.noweight.RData'))
        #get.stats(N=N,m=m,alpha=alpha,h2=h2,times.all=20,inf=inf)
        #}
      }
    }
  }
}















#apply(res,2,mean)
apply(res,2,sd)


path0 <- '/home/songs/genomic_control/20201218simu/new/revision220127_fixR/'
setwd(path0)
for(h2 in c(0.2,0.5,0.05,0.1)){
  for(alpha in c(0.005,0.01,0.05,0.1)){
    for(m in c(1e5)){
      #  for(N in c(2000,5000,10000,2e4,5e4,1e5,2e5)){
      for (inf in c(1,1.1)){
        load(paste0(path0,'/plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
                    '_inf_',inf,'.plot_dat50.new_2stage_inf1stage.RData'))
        ## hess
        
        # res.all <- c()
        h2.est <- a.est <- c()
        N <- rep(c(5000,1e4,2e4,5e4),each=50)
        for(nn in c(5000,1e4,2e4,5e4)){
          path <- paste0(path0,'/hess/sumstats/inf_',inf,'/h2_',h2,'/alpha_',alpha,'/m_',m,'/N_',nn)
          load(paste0(path,'/res.RData'))
          h2.est <- c(h2.est,unlist(lapply(res,'[',2)))
          #res.all <- get.a(n.gwas,m,alpha=alpha,inf=inf)  
          a.est <- c(a.est,unlist(lapply(res,'[',1)))
        }
        Method <- rep('HESS',50*4)
        dat <- rbind(dat,data.frame(Method,N,h2.est,a.est))
        
        ### hdl
        h20 <- h2
        h2.est <- a.est <- c()
        N <- rep(c(5000,1e4,2e4,5e4),each=50)
        for(nn in c(5000,1e4,2e4,5e4)){
          path <- paste0(path0,'/inf_',inf,'/h2_',h20,'/alpha_',alpha,'/m_',m,'/N_',nn)
          load(paste0(path,'/hdl_res.RData'))
          h2.est <- c(h2.est,h2)
          h2 <- h20
          #res.all <- get.a(n.gwas,m,alpha=alpha,inf=inf)  
          a.est <- c(a.est,rep(1,50))
        }
        Method <- rep('HDL',50*4)
        dat <- rbind(dat,data.frame(Method,N,h2.est,a.est))
        
        dir.create(paste0(path0,'/plot_dat'))
        save(dat,file=paste0(path0,'/plot_dat/h2_',h20,'_alpha_',alpha,'_m_',m,
                             '_inf_',inf,'.plot_dat50.new_2stage_inf1stage.all.RData'))
        #get.stats(N=N,m=m,alpha=alpha,h2=h2,times.all=20,inf=inf)
        #}
      }
    }
  }
}



#### plot

library(ggsci)
library(ggplot2)
#'{}^0*italic(D)~plain((rarefied))'

get.plot.both <- function(h2,alpha,m,inf){
  library(scales)
  col <- pal_npg("nrc")(2)
  load(paste0('./plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
              '_inf_',inf,'.plot_dat50.new_2stage_inf1stage.all.RData'))
  # dat <- dat[dat$N!=1e5,] ## not include 1e5 here
  # dat <- dat[dat$N!=2000,] ## not include 1e5 here
  dat1 <- dat[dat$Method=='LDER',]
  
  dim(dat1)
  dat1$Method <- 'LDER (optimal weights)'
  load(paste0('./plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
              '_inf_',inf,'.plot_dat50.noweight.RData'))
  dat2 <- dat[dat$Method=='LDER',]
  dim(dat2)
  dat2$Method <- 'LDER (no weighting)'
  dat <- rbind(dat1,dat2)
  
  dat$N1 <- paste0('n = ',dat$N)
  dat$N1 <- factor(dat$N1,levels=unique(dat$N1))
  
   p1 <-  ggplot(dat, aes(x = Method, y = h2.est)) + 
    geom_boxplot(aes(fill = Method),
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = col)+
    #scale_fill_npg()+scale_color_npg()+
    theme_bw()+
    #theme(legend.position = 'none')+ 
    xlab(' ')+
    ylab(expression(Estimated~h^2))+
    theme(panel.grid=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line=element_blank(),
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
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
   # if(inf==1){
  # p1 <- p1+labs(title=substitute(list(h^2 == h2, alpha == aa, 'no inflation'), list(h2=h2, aa=alpha)))+
  #   theme(plot.title = element_text(hjust = 0.5))
  # }else{
  #   p1 <- p1+labs(title=substitute(list(h^2 == h2, alpha == aa, 'Inf'==inf), list(h2=h2, aa=alpha, inf=inf)))+
  #     theme(plot.title = element_text(hjust = 0.5))
  #   
  # }
  p2 <- ggplot(dat, aes(x = Method, y = a.est)) + 
    geom_boxplot(aes(fill = Method),
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = col)+
    #scale_fill_npg()+scale_color_npg()+
    theme_bw()+
    #theme(legend.position = 'none')+ 
    xlab(' ')+
    ylab(expression(Estimated~inflation))+
    theme(panel.grid=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line=element_blank(),
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
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
  p2
  if(inf>1){
    title <- grid::textGrob(substitute(list(h^2 == h2, alpha == aa, 'Inf'==inf), list(h2=h2, aa=alpha, inf=inf)), gp=grid::gpar(fontsize=13, fontface='bold'))
  }else{
    title <- grid::textGrob(substitute(list(h^2 == h2, alpha == aa, 'no inflation'), list(h2=h2, aa=alpha)), gp=grid::gpar(fontsize=13, fontface='bold'))
  }
  p3 <- lemon::grid_arrange_shared_legend(p1,p2,nrow=2,ncol=1,top=title)
  p3
  return(p3)
}
## two h2
get.plot.h2 <- function(h2,alpha,m,inf){
  library(scales)
  col <- pal_npg("nrc")(2)
  load(paste0('./plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
              '_inf_',inf,'.plot_dat50.new_2stage_inf1stage.all.RData'))
  # dat <- dat[dat$N!=1e5,] ## not include 1e5 here
  # dat <- dat[dat$N!=2000,] ## not include 1e5 here
  dat1 <- dat[dat$Method=='LDER',]
  
  dim(dat1)
  dat1$Method <- 'LDER (optimal weights)'
  load(paste0('./plot_dat/h2_',h2,'_alpha_',alpha,'_m_',m,
              '_inf_',inf,'.plot_dat50.noweight.RData'))
  dat2 <- dat[dat$Method=='LDER',]
  dim(dat2)
  dat2$Method <- 'LDER (no weighting)'
  dat <- rbind(dat1,dat2)
  
  dat$N1 <- paste0('n = ',dat$N)
  dat$N1 <- factor(dat$N1,levels=unique(dat$N1))
  dat$Method <- factor(dat$Method,levels=unique(dat$Method))
  p1 <-  ggplot(dat, aes(x = Method, y = h2.est)) + 
    geom_boxplot(aes(fill = Method),
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = col)+
    #scale_fill_npg()+scale_color_npg()+
    theme_bw()+
    #theme(legend.position = 'none')+ 
    xlab(' ')+
    ylab(expression(Estimated~h^2))+
    theme(panel.grid=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line=element_blank(),
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
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
  if(inf>1){
    title <- grid::textGrob(substitute(list(h^2 == h2, alpha == aa, 'Inf'==inf), list(h2=h2, aa=alpha, inf=inf)), gp=grid::gpar(fontsize=13, fontface='bold'))
  }else{
    title <- grid::textGrob(substitute(list(h^2 == h2, alpha == aa, 'no inflation'), list(h2=h2, aa=alpha)), gp=grid::gpar(fontsize=13, fontface='bold'))
  }
  p3 <- lemon::grid_arrange_shared_legend(p1,nrow=1,ncol=1,top=title)
  p3
  return(p3)
}

# title <- grid::textGrob(substitute(list(h^2 == h2, alpha == aa, 'Inf'==inf), list(h2=h2, aa=alpha, inf=inf)), gp=grid::gpar(fontsize=14, fontface='bold'))
# p3 <- lemon::grid_arrange_shared_legend(p1,p2,nrow=2,ncol=1,
#                                         top=title)
# # p1 + labs(title=bquote(alpha*(.(h2)) ","*nu))+
# #   theme(plot.title = element_text(hjust = 0.5))
# p3+labs(title=substitute(list(h^2 == h2, alpha == aa, 'Inf'==inf), list(h2=h2, aa=alpha, inf=inf)))+
# theme(plot.title = element_text(hjust = 0.5))
# # p1 <- get.plot(h2=0.5,alpha=0.01,m=1e5,inf=1,type='h2')
# # p1
# # p2 <- get.plot(h2=0.5,alpha=0.01,m=1e5,inf=1,type='inf')
# # p2

p3 <- get.plot.h2(h2=0.5,alpha=0.05,m=1e5,inf=1)
p3
pdf('no_weight.pdf',height=3.3,width = 9)
get.plot.h2(h2=0.5,alpha=0.05,m=1e5,inf=1)
get.plot.h2(h2=0.5,alpha=0.05,m=1e5,inf=1.1)
dev.off()


pdf('wider_range_inf1.1.pdf',height=5.5,width = 9)
get.plot.both(h2=0.05,alpha=0.05,m=1e5,inf=1)
get.plot.both(h2=0.1,alpha=0.05,m=1e5,inf=1)
get.plot.both(h2=0.05,alpha=0.1,m=1e5,inf=1)
get.plot.both(h2=0.1,alpha=0.1,m=1e5,inf=1)
dev.off()



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









