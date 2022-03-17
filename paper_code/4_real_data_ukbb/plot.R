
#######
path <- '/home/songs/UKB_summary/20210211/'
setwd(path)

library(data.table)
dat <- fread('0224dat.analyse.txt',sep='\t',head=T)
dat <- dat[dat$pheno.code!='20122',]
## File: file name
## 
dat <- as.data.frame(dat)
dat <- dat[!duplicated(dat$pheno.code),]
dat.quan <- dat[dat$variable_type!='binary',]
dat.qual <- dat[dat$variable_type=='binary',]
#dat <- dat.quan
h2.ldsc <- h2.lder <- h2.ldsc.sd <- h2.lder.sd <- inf.lder <- 
  inf.ldsc <- inf.lder.sd <- inf.ldsc.sd <- c()
for(w in 1:nrow(dat)){
  setwd(path)
  print(w)
  pheno <- dat$pheno.code[w]
  if(dat$variable_type[w]=='binary'){
    setwd('/home/songs/UKB_summary/20210211/res_ukb/newinf_qualitative')
    load(paste0(pheno,'_res_220212.RData'))
    h2.lder <- c(h2.lder,res$lder$h2.lia)
    h2.lder.sd <- c(h2.lder.sd,res$lder$h2.sd.lia)
    inf.lder <- c(inf.lder,res$lder$inf)
    inf.lder.sd <- c(inf.lder.sd,res$lder$inf.sd)
    
    h2.ldsc <- c(h2.ldsc,res$ldsc$h2.lia)
    h2.ldsc.sd <- c(h2.ldsc.sd,res$ldsc$h2.sd.lia)
    inf.ldsc <- c(inf.ldsc,res$ldsc$inf)
    inf.ldsc.sd <- c(inf.ldsc.sd,res$ldsc$inf.sd)
    
  }else{
    
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
  
  
}


dat1 <- data.frame(pheno=dat$pheno.code,
                   lder=h2.lder,
                   lder.sd=h2.lder.sd,
                   ldsc=h2.ldsc,
                   ldsc.sd=h2.ldsc.sd,
                   variable_type=dat$variable_type,
                   category=dat$`Main Categoriy`,
                   description=dat$description)
dat1$pv.lder <- pnorm(-(dat1$lder)/dat1$lder.sd)
dat1$pv.ldsc <- pnorm(-(dat1$ldsc)/dat1$ldsc.sd)
dat1$pv.diff <- pnorm(-abs(dat1$ldsc-dat1$lder)/
                        sqrt(dat1$lder.sd^2+dat1$ldsc.sd^2))

length(which(dat1$pv.diff<0.05/814))
length(which(dat1$pv.lder<0.05/814))
length(which((dat1$pv.lder<0.05/814)&(dat1$pv.ldsc>=0.05/814)))
setwd('/home/songs/UKB_summary/20210211/res_ukb/')
write.table(dat1,'result_type_newinf.txt',quote=F,row.names=F,col.names=T,sep='\t')








##### plot




library(data.table)
dat <-   fread('results/result_type_newinf.txt',sep='\t')
head(dat)  
library(ggplot2)
dat <- dat[which(dat$pheno!='20122'),]
get.figure <- function(dat,binary=T){
  label <- gsub(".*: ","",dat$description)
  label[which(label=='None of the above')] <- dat$description[which(label=='None of the above')]
  dat$description <- label
  if(binary){
    dat1 <- dat[dat$variable_type=='binary',]
    title <- 'Dichotomous Phenotypes'
  }else{
    dat1 <- dat[dat$variable_type!='binary',]
    title <- 'Quantitative Phenotypes'
  }
  dat2 <- dat1[dat1$pv.diff<0.05/nrow(dat),]
  dat3 <- dat1[which((dat1$pv.lder<0.05/nrow(dat))&(dat1$pv.ldsc>=0.05/nrow(dat))),]
  
  
  dat2 <- dat2[order(dat2$pv.lder),]
  dat4 <- dat2[1:10,]
  p1 <- ggplot(dat1, aes(x = ldsc, y = lder)) +
    geom_point(size =4,aes(color=pv.diff)) +
    scale_color_continuous(trans = 'reverse')+
    geom_point(aes(x=ldsc,y=lder),size =4.5,shape=1,colour='orange',dat2) +
    geom_point(aes(x=ldsc,y=lder),size =4.5,shape=1,colour='#d42020',dat3) +
    geom_point(aes(x=ldsc,y=lder),size =5,shape=1,colour='#d42020',dat3) +
    geom_point(aes(x=ldsc,y=lder),size =5,shape=1,colour='orange',dat2) +
    #geom_errorbar(aes(ymin = lder-lder.sd, ymax =  lder+lder.sd,color=-pv.lder),dat2)+
    # geom_errorbar(aes(xmin = ldsc-ldsc.sd, xmax =  ldsc+ldsc.sd,color=-pv.ldsc),dat2)+
    theme(
      panel.background = element_rect(fill = "white",colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size=80,color="black",family="sans"),
      axis.title = element_text(size=90,face = "bold",family="sans"),
      legend.text =  element_text(size=80,family="sans"),
      legend.title = element_text(size=90,family="sans"),
      legend.position='none',
      # axis.title.x = element_text(hjust=5),
      axis.ticks = element_blank())+
    geom_abline( slope=1,
                 intercept=0,colour='#E41A1C',lty=2)+
    xlab(expression(h^2~'estimated by LDSC'))+ylab(expression(h^2~'estimated by LDER'))
  library(ggrepel)
  p1
  p2 <- p1+  geom_label_repel(data=dat4,aes(label = description),
                              box.padding   = 0.4,
                              point.padding = 0.5,
                              label.size = 0.1,
                              #ylim=c(0.01,NA),
                              direction="both",
                              nudge_y = 0.01,
                              #nudge_x = 0.15,
                              segment.color = 'grey50',force=1) +
    theme_classic()+
    labs(color = "p value")+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5))
  p2
  return(p2)
}
get.figure.new1 <- function(dat,binary=T){
  label <- gsub(".*: ","",dat$description)
  label[which(label=='None of the above')] <- dat$description[which(label=='None of the above')]
  label1 <- unlist(lapply(strsplit(label,' /'),'[',1))
  dat$description <- label1
  if(binary){
    dat1 <- dat[dat$variable_type=='binary',]
    title <- 'Dichotomous Phenotypes'
  }else{
    dat1 <- dat[dat$variable_type!='binary',]
    title <- 'Quantitative Phenotypes'
  }
  dat2 <- dat1[dat1$pv.diff<0.05/nrow(dat),]
  dat3 <- dat1[which((dat1$pv.lder<0.05/nrow(dat))&(dat1$pv.ldsc>=0.05/nrow(dat))),]
  
  
  dat2 <- dat2[order(dat2$pv.lder),]
  rr <- min(nrow(dat2),6)
  dat4 <- dat2[1:rr,]
  p1 <- ggplot(dat1, aes(x = ldsc, y = lder)) +
    geom_point(size =4,aes(color=pv.diff)) +
    scale_color_continuous(trans = 'reverse')+
    geom_point(aes(x=ldsc,y=lder),size =4.5,shape=1,colour='orange',dat2) +
    # geom_point(aes(x=ldsc,y=lder),size =4.5,shape=1,colour='#d42020',dat3) +
    # geom_point(aes(x=ldsc,y=lder),size =5,shape=1,colour='#d42020',dat3) +
    # geom_point(aes(x=ldsc,y=lder),size =5,shape=1,colour='orange',dat2) +
    #geom_errorbar(aes(ymin = lder-lder.sd, ymax =  lder+lder.sd,color=-pv.lder),dat2)+
    # geom_errorbar(aes(xmin = ldsc-ldsc.sd, xmax =  ldsc+ldsc.sd,color=-pv.ldsc),dat2)+
    theme(
      panel.background = element_rect(fill = "white",colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size=8,color="black",family="sans"),
      axis.title = element_text(size=9,face = "bold",family="sans"),
      legend.text =  element_text(size=8,family="sans"),
      legend.title = element_text(size=9,family="sans"),
      legend.position='none',
      # axis.title.x = element_text(hjust=5),
      axis.ticks = element_blank())+
    geom_abline( slope=1,
                 intercept=0,colour='#E41A1C',lty=2)+
    xlab(expression(h^2~'estimated by LDSC'))+ylab(expression(h^2~'estimated by LDER'))
  library(ggrepel)
  p1
  p2 <- p1+  geom_label_repel(data=dat4,aes(label = description),
                              box.padding   = 0.4,
                              point.padding = 0.5,
                              label.size = 0.1,
                              #ylim=c(0.01,NA),
                              direction="both",
                              nudge_y = 0.01,
                              #nudge_x = 0.15,
                              max.overlaps=100,
                              segment.color = 'grey50',force=1) +  
    theme_classic()+
    labs(color = "p value")+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5))
  p2
  return(p2)
}
get.figure.new1(dat,binary=T)
get.figure.new2 <- function(dat,binary=T){
  label <- gsub(".*: ","",dat$description)
  label[which(label=='None of the above')] <- dat$description[which(label=='None of the above')]
  dat$description <- label
  if(binary){
    dat1 <- dat[dat$variable_type=='binary',]
    title <- 'Dichotomous Phenotypes'
  }else{
    dat1 <- dat[dat$variable_type!='binary',]
    title <- 'Quantitative Phenotypes'
  }
  dat2 <- dat1[dat1$pv.diff<0.05/nrow(dat),]
  dat3 <- dat1[which((dat1$pv.lder<0.05/nrow(dat))&(dat1$pv.ldsc>=0.05/nrow(dat))),]
  
  
  dat2 <- dat2[order(dat2$pv.lder),]
  dat3 <- dat3[order(dat3$pv.lder),]
  dat4 <- dat3[1:6,]
  p1 <- ggplot(dat3, aes(x = ldsc, y = lder))+
    #   geom_point(size =2,color='gray') +
    # p1 <- ggplot(dat1, aes(x = ldsc, y = lder)) +
    geom_point(aes(x=ldsc,y=lder),size =3,colour='gray',dat1) +
    geom_point(size =3,aes(color=pv.lder)) +
    #scale_color_continuous(trans = 'reverse')+
    scale_color_gradient2(low="#ffd5d5",mid='#ff5454', high="#bc0303")+
    #geom_point(aes(x=ldsc,y=lder),size =2.5,shape=1,colour='orange',dat2) +
    
    # geom_point(aes(x=ldsc,y=lder),size =5,shape=1,colour='#d42020',dat3) +
    # geom_point(aes(x=ldsc,y=lder),size =5,shape=1,colour='orange',dat2) +
    #geom_errorbar(aes(ymin = lder-lder.sd, ymax =  lder+lder.sd,color=-pv.lder),dat2)+
    # geom_errorbar(aes(xmin = ldsc-ldsc.sd, xmax =  ldsc+ldsc.sd,color=-pv.ldsc),dat2)+
    theme(
      panel.background = element_rect(fill = "white",colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size=8,color="black",family="sans"),
      axis.title = element_text(size=9,face = "bold",family="sans"),
      legend.text =  element_text(size=8,family="sans"),
      legend.title = element_text(size=9,family="sans"),
      legend.position='none',
      # axis.title.x = element_text(hjust=5),
      axis.ticks = element_blank())+
    geom_abline( slope=1,
                 intercept=0,colour='#E41A1C',lty=2)+
    xlab(expression(h^2~'estimated by LDSC'))+ylab(expression(h^2~'estimated by LDER'))
  library(ggrepel)
  p1
  p2 <- p1+  geom_label_repel(data=dat4,aes(label = description),
                              box.padding   = 0.4,
                              point.padding = 0.5,
                              label.size = 0.1,
                              #ylim=c(0.01,NA),
                              direction="both",
                              nudge_y = 0.05,
                              max.overlaps=100,
                              #nudge_x = 0.15,
                              segment.color = 'grey50',force=1) +
    theme_classic()+
    labs(color = "p value")+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5))
  p2
  return(p2)
}
get.figure.new2(dat,binary=T)
library(ggplot2)
p1 <- get.figure.new1(dat,binary=F)
p1
p2 <- get.figure.new1(dat,binary=T)
p2
p3 <- get.figure.new2(dat,binary=F)
p4 <- get.figure.new2(dat,binary=T)
pdf('./ukb_alltraits/ukb_traits_v61.pdf',height=6,width = 13.5,onefile = F)
lemon::grid_arrange_shared_legend(p1,p2,nrow=1,ncol=2)
dev.off()
pdf('./ukb_alltraits/ukb_traits_v62.pdf',height=6,width = 13.5,onefile = F)
lemon::grid_arrange_shared_legend(p3,p4,nrow=1,ncol=2)
dev.off()




dat[dat$pv.diff<0.05/814,]




library(data.table)
dat <-   fread('results/result_type_newinf.txt',sep='\t')
head(dat)  
library(ggplot2)
dat <- dat[which(dat$pheno!='20122'),]
dim(dat)
dat2 <- dat[dat$pv.diff<0.05/814,]
dim(dat2)
dat3 <- dat2[,c('pheno','description','variable_type')]
dat3$lder <- paste0(format(round(dat2$lder,3), nsmall = 3),' (',format(round(dat2$lder.sd,3), nsmall = 3),')')
dat3$ldsc <- paste0(format(round(dat2$ldsc,3), nsmall = 3),' (',format(round(dat2$ldsc.sd,3), nsmall = 3),')')
head(dat3)  
dat3$pv <- dat2$pv.diff
write.csv(dat3,'table2.csv',quote=T,row.names=F)

dat4 <- dat[(dat$pv.lder<0.05/814)&(dat$pv.ldsc>=0.05/814),]
dim(dat4)
dat4[1,]
h2.lder <- dat4$lder
h2.lder.sd <- dat4$lder.sd
h2.ldsc <- dat4$ldsc
h2.ldsc.sd <- dat4$ldsc.sd
dat4$lder.h2 <- paste0(format(round(h2.lder,3), nsmall = 3),' (',format(round(h2.lder.sd,3), nsmall = 3),')')
dat4$ldsc.h2 <- paste0(format(round(h2.ldsc,3), nsmall = 3),' (',format(round(h2.ldsc.sd,3), nsmall = 3),')')


write.csv(dat4,'table5.csv',quote=T,row.names=F)


