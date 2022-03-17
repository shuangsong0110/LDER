library(data.table)
library(scales)
squish_trans <- function(from, to, factor) {
  
  trans <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squished", trans, inv))
}

dat <- fread('./results/res_10traits_newinf.csv')
head(dat)
type <- c('Quantitative Phenotypes','Dichotomous Phenotypes')
library(ggplot2)
library(ggsci)
Trait <- rep(c('LDL','HDL','TG','TC','HGT','BMI'),5)
Trait <- factor(Trait,levels=unique(Trait))
Method <- rep(c('BOLT-LMM','LDER','LDSC','HESS','HDL'),each=6)
h2 <- c(t(dat[1:6,2]),t(dat[1:6,4]),t(dat[1:6,6]),t(dat[1:6,8]),t(dat[1:6,10]))
sd <- c(t(dat[1:6,3]),t(dat[1:6,5]),t(dat[1:6,7]),t(dat[1:6,9]),t(dat[1:6,11]))
dat1 <- data.frame(Trait,Method,h2,sd)
head(dat1)
dat1$sd1 <- dat1$h2-dat1$sd
dat1$sd2 <- dat1$h2+dat1$sd
col <- pal_npg("nrc")(5)[c(4,1,2,3,5)]
#col <- c('#ffdb95',pal_npg("nrc")(4))
library(ggsci)
dat1$Method <- factor(dat1$Method,levels=unique(dat1$Method))
p1 = ggplot(data = dat1,
            aes(x =Trait, y = h2, fill = Method))+
  geom_bar(position = 'dodge',stat = 'identity',width = 0.85)+
  xlab(expression(Trait))+
  #xlab(expression(p~value~bins~'in'~discovery~cohort~'('-log[10]~')'))+
  ylab(expression(Estimated~h^2))+
  geom_errorbar(aes(x=Trait,y=h2,ymin=sd1, ymax=sd2),width=.3,size=.8,
                position=position_dodge(.85),colour=rgb(136/255,136/255,136/255))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text=element_text(size=14),
        axis.line=element_line(size=1),
        legend.title =  element_text(size=11),
        legend.text =  element_text(size=12),
        axis.title=element_text(size=16))+
  theme(legend.position = c(0.7,0.7),
        legend.title = element_blank(),
        #  legend.title = element_text(size=16,face="bold"),
        legend.text = element_text(size=16))+
  ggtitle(type[1])+theme(plot.title = element_text(hjust = 0.45,size = 17, face = "bold"))+
#  scale_fill_npg()
  scale_fill_manual(values = col)

p1

dat <- fread('./results/res_10traits_newinf.csv')
Trait <- rep(c('ATH','CAD','SCZ','T2D'),5)
Trait <- factor(Trait,levels=unique(Trait))
Method <- rep(c('BOLT-LMM','LDER','LDSC','HESS','HDL'),each=4)
h2 <- c(t(dat[7:10,2]),t(dat[7:10,4]),t(dat[7:10,6]),t(dat[7:10,8]),t(dat[7:10,10]))
sd <- c(t(dat[7:10,3]),t(dat[7:10,5]),t(dat[7:10,7]),t(dat[7:10,9]),t(dat[7:10,11]))
dat1 <- data.frame(Trait,Method,h2,sd)
head(dat1)
dat1$sd1 <- dat1$h2-dat1$sd
dat1$sd2 <- dat1$h2+dat1$sd
dat1$sd1[dat1$sd1<0] <- 0
library(ggsci)
library(scales)
dat1$Method <- factor(dat1$Method,levels=unique(dat1$Method))
p2 = ggplot(data = dat1,
            aes(x =Trait, y = h2, fill = Method))+
  geom_bar(position = 'dodge',stat = 'identity',width = .75)+
  xlab(expression(Trait))+
  #xlab(expression(p~value~bins~'in'~discovery~cohort~'('-log[10]~')'))+
  ylab(expression(Estimated~h^2))+
  geom_errorbar(aes(x=Trait,y=h2,ymin=sd1, ymax=sd2),width=.3,size=.8,
                position=position_dodge(.75),colour=rgb(136/255,136/255,136/255))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text=element_text(size=14),
        axis.line=element_line(size=1),
        legend.title =  element_text(size=11),
        legend.text =  element_text(size=12),
        axis.title=element_text(size=16))+
  theme(legend.position = c(0.7,0.7),
        legend.title = element_blank(),
        #  legend.title = element_text(size=16,face="bold"),
        legend.text = element_text(size=16))+
  ggtitle(type[2])+theme(plot.title = element_text(hjust = 0.45,size = 17, face = "bold"))+
  scale_fill_manual(values = col)+
  #coord_cartesian(ylim = c(0, 1))+
  scale_y_continuous(trans = squish_trans(0.3,1.8, 40),
                     breaks = c(0, 0.1, 0.2,0.3,1.8,1.9,2))

p2

lemon::grid_arrange_shared_legend(p1,p2)



# p2
dir.create('./plot')
pdf('./plot/ukb_10traits_all7.pdf',height = 5,width = 12,onefile = F)
lemon::grid_arrange_shared_legend(p1,p2)
dev.off()
## old version
if(F){
############################
#trait <- c('BC','CD','PrCa','RA','UC')
library(data.table)
dat <- fread('./results/res_10traits.csv')
head(dat)
type <- c('Quantitative Phenotypes','Dichotomous Phenotypes')
library(ggplot2)
library(ggsci)
Trait <- rep(c('LDL','HDL','TG','TC','HGT','BMI'),3)
Trait <- factor(Trait,levels=unique(Trait))
Method <- rep(c('BOLT-LMM','LDER','LDSC'),each=6)
h2 <- c(t(dat[1:6,2]),t(dat[1:6,4]),t(dat[1:6,6]))
sd <- c(t(dat[1:6,3]),t(dat[1:6,5]),t(dat[1:6,7]))
dat1 <- data.frame(Trait,Method,h2,sd)
head(dat1)
dat1$sd1 <- dat1$h2-dat1$sd
dat1$sd2 <- dat1$h2+dat1$sd
library(ggsci)
p1 = ggplot(data = dat1,
            aes(x =Trait, y = h2, fill = Method))+
  geom_bar(position = 'dodge',stat = 'identity',width = 0.6)+
  xlab(expression(Trait))+
  #xlab(expression(p~value~bins~'in'~discovery~cohort~'('-log[10]~')'))+
  ylab(expression(Estimated~h^2))+
  geom_errorbar(aes(x=Trait,y=h2,ymin=sd1, ymax=sd2),width=.3,size=.8,
                position=position_dodge(.6),colour=rgb(136/255,136/255,136/255))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text=element_text(size=14),
        axis.line=element_line(size=1),
        legend.title =  element_text(size=11),
        legend.text =  element_text(size=12),
        axis.title=element_text(size=16))+
  theme(legend.position = c(0.7,0.7),
        legend.title = element_blank(),
        #  legend.title = element_text(size=16,face="bold"),
        legend.text = element_text(size=16))+
  ggtitle(type[1])+theme(plot.title = element_text(hjust = 0.45,size = 17, face = "bold"))+scale_fill_npg()


p1

dat <- fread('./results/res_10traits.csv')
Trait <- rep(c('ATH','CAD','SCZ','T2D'),3)
Trait <- factor(Trait,levels=unique(Trait))
Method <- rep(c('BOLT-LMM','LDER','LDSC'),each=4)
h2 <- c(t(dat[7:10,2]),t(dat[7:10,4]),t(dat[7:10,6]))
sd <- c(t(dat[7:10,3]),t(dat[7:10,5]),t(dat[7:10,7]))
dat1 <- data.frame(Trait,Method,h2,sd)
head(dat1)
dat1$sd1 <- dat1$h2-dat1$sd
dat1$sd2 <- dat1$h2+dat1$sd
dat1$sd1[dat1$sd1<0] <- 0
library(ggsci)
p2 = ggplot(data = dat1,
            aes(x =Trait, y = h2, fill = Method))+
  geom_bar(position = 'dodge',stat = 'identity',width = 0.6)+
  xlab(expression(Trait))+
  #xlab(expression(p~value~bins~'in'~discovery~cohort~'('-log[10]~')'))+
  ylab(expression(Estimated~h^2))+
  geom_errorbar(aes(x=Trait,y=h2,ymin=sd1, ymax=sd2),width=.3,size=.8,
                position=position_dodge(.6),colour=rgb(136/255,136/255,136/255))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text=element_text(size=14),
        axis.line=element_line(size=1),
        legend.title =  element_text(size=11),
        legend.text =  element_text(size=12),
        axis.title=element_text(size=16))+
  theme(legend.position = c(0.7,0.7),
        legend.title = element_blank(),
        #  legend.title = element_text(size=16,face="bold"),
        legend.text = element_text(size=16))+
  ggtitle(type[2])+theme(plot.title = element_text(hjust = 0.45,size = 17, face = "bold"))+scale_fill_npg()

p2

lemon::grid_arrange_shared_legend(p1,p2)



# p2
pdf('./plot/ukb_10traits.pdf',height = 5,width = 10,onefile = F)
lemon::grid_arrange_shared_legend(p1,p2)
dev.off()
}
