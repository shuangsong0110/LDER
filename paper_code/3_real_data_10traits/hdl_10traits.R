library(HDL)


data(gwas1.example)
LD.path <- '/home/songs/UKB_data_qc2/ukb_hdl/'
LD.path <- '/home/songs/UKB_data_qc2/ukb_hdl_imp/UKB_imputed_SVD_eigen99_extraction/'
#LD.path <- "/Path/to/reference/UKB_array_SVD_eigen90_extraction"
# res.HDL <- HDL.h2(gwas.df = gwas1.example, LD.path = LD.path)
# res.HDL




library(data.table)
for(trait in c('T2D','ASTHMA','CAD','SCZ', 'HEIGHT','BMI', 'LDL','TG','TC','HDL')){
  load(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait,'/lder.assoc.RData'))
  path <- paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait,'/HDL/')
  dir.create(path,recursive = T)
  colnames(assoc) <- c('BP','SNP','A1','A2','CHR','Z','N')
  res <- HDL.h2(gwas.df = assoc, LD.path = LD.path)
  save(res,file=paste0(path,'/res.RData'))
  #write.table(dat,paste0(path,'/summs.txt'),quote=F,row.names=F,col.names=T)
  print(trait)
}


for(trait in c( 'HEIGHT','BMI', 'LDL','TG','TC','HDL')){
  #load(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait,'/lder.assoc.RData'))
  path <- paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait,'/HDL/')
 load(paste0(path,'/res.RData'))
  #write.table(dat,paste0(path,'/summs.txt'),quote=F,row.names=F,col.names=T)
  print(trait)
  print(res$h2)
  print(res$h2.se)
}

cvt.h2 <- function(K,P,h2){
  # K: prevalence
  # P: proportion of cases
  zv <- dnorm(qnorm(K))
  h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
  return(h2_liab)
}


for(trait in c( 'T2D','ASTHMA','CAD','SCZ')){
  #load(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait,'/lder.assoc.RData'))
  path <- paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait,'/HDL/')
  load(paste0(path,'/res.RData'))
  #write.table(dat,paste0(path,'/summs.txt'),quote=F,row.names=F,col.names=T)
  print(trait)
  setwd(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait))
  phen <- fread('plink.log',fill=T)
  ncase <- as.numeric(phen$V4[25])
  ncontrol <- as.numeric(phen$V4[24])-ncase
  K <- P <- ncase/(ncase+ncontrol)
  h2.lia <- cvt.h2(K,P,res$h2)
  h2.sd.lia <- cvt.h2(K,P,res$h2.se)
  
  print(h2.lia)
  print(h2.sd.lia)
}



















if(F){


setwd('/home/songs/genomic_control/hdl/1kg_wtccc')
### 
#wtccc_1kg overlap SNPs::
chr <- 22
library(data.table)
bim1 <- fread(paste0('/home/songs/summerResearch/main/reference_1KG/1000G.EUR.QC.',chr,'.bim'))
bim2 <- fread(paste0('/home/songs/CV3/newdata_2020/WTCCC/RA/chr/chr',chr,'.bim'))
bim1$order <- 1:nrow(bim1)
bim <- merge(bim1,bim2,by='V2')
bim <- bim[order(bim[,4]),]
order <- bim$order
bim <- bim[,c(2,1,3:6)]

snps.name.list <- bim$V2
save(snps.name.list,file='UKB_snp_list_array.vector_form.RData')
nsnps.list <- list()
for(k in 1:21){nsnps.list[[k]] <- NULL}
nsnps.list[[22]] <- nrow(bim)
save(nsnps.list,file='UKB_snp_counter_array.RData')
write.table(bim,'ukb_chr22.1_n336000.bim',row.names = F,col.names = F,quote=F)

library(plink2R)
dat <- read_plink(paste0('/home/songs/summerResearch/main/reference_1KG/1000G.EUR.QC.',chr))
n <- nrow(dat$fam)
X <- dat$bed[,order]
R <- cor(X)
temp <- eigen(R)

source("/home/songs/genomic_control/Rcode/ldsr_gai_new.R")
LDsc <- ldscore(R,n)
V <- temp$vectors
lam <- temp$values

save(LDsc,V,lam,file='ukb_array_chr22.1_n336000_500banded_90eigen.rda')


}



























