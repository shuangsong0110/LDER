### write summary statistics
library(data.table)
for(trait in c('T2D','ASTHMA','CAD','SCZ', 'HEIGHT','BMI', 'LDL','TG','TC','HDL')){
  load(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait,'/lder.assoc.RData'))
  path <- paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait,'/hess/')
  dir.create(path,recursive = T)
  dat <- assoc[,c(2,5,1,3,4,6,7)]
  colnames(dat) <- c('SNP','CHR','BP','A1','A2','Z','N')
  write.table(dat,paste0(path,'/summs.txt'),quote=F,row.names=F,col.names=T)
  print(trait)
}



### UKB ld
cd /home/songs/genomic_control/hess/hess-0.5.3-beta/
path0='/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/'
for trait in 'T2D' 'ASTHMA' 'CAD' 'SCZ' 'HEIGHT' 'BMI' 'LDL' 'TG' 'TC' 'HDL'
do
mkdir ${path0}${trait}/hess/result
for chrom in {1..22}
do
python2 ./hess.py \
--local-hsqg ${path0}${trait}/hess/summs.txt \
--chrom ${chrom} \
--bfile /home/songs/genomic_control/ukb/dat/alldat_sam1e4.qc/chr/chr${chrom} \
--partition /home/songs/jun/ldetect-data/UKB_hess/EUR/fourier_ls-chr${chrom}.bed \
--out ${path0}${trait}/hess/result/step1
done
python2 ./hess.py --prefix ${path0}${trait}/hess/result/step1 --out ${path0}${trait}/hess/result/step2 --reinflate-lambda-gc 1.1 --path ${path0}${trait}/hess/result/
done
done


## lambda gc
cd /home/songs/genomic_control/hess/hess-0.5.3-beta/
  path0='/home/songs/genomic_control/ukb/analysis/ukb_ld/hess/'
for N in '5000' '10000' '20000' '50000'
do
h2=0.5
alpha=0.01
for times in {1..50}
do
python2 ./hess.py --prefix ${path0}/result_ukb/h2_${h2}/alpha_${alpha}/N_${N}/step1 --out ${path0}/result_ukb/h2_${h2}/alpha_${alpha}/N_${N}/step2_times${times} --reinflate-lambda-gc 1.1 --path ${path0}/result_ukb/h2_${h2}/alpha_${alpha}/N_${N}/times${times}_
done
done


cvt.h2 <- function(K,P,h2){
  # K: prevalence
  # P: proportion of cases
  zv <- dnorm(qnorm(K))
  h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
  return(h2_liab)
}


for(trait in c('T2D','ASTHMA','CAD','SCZ')){
  setwd('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/')
  dir.create(trait)
  print(trait)
  setwd(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait))
 # load('lder.assoc.RData')
  ldpath <- '/home/songs/genomic_control/ukb/20210204/ukb/LD/'
  ldpath.shrink <- '/home/songs/genomic_control/ukb/20210204/ukb/LD.shrink/'
  
  # res <- ld.all(assoc,assoc$N[1],ldpath=ldpath,ldpath.shrink=ldpath.shrink,method='both',
  #               type='jack',rough=F,twostage=T,cores=40)
  res <- fread(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',
                      trait,'/hess/result/step2.log'))
  res <- as.character(res[17,1])
  h2 <- as.numeric(strsplit(res,' ')[[1]][5])
  sd <- strsplit(res,' ')[[1]][6]
  sd <- strsplit(sd,')')[[1]][1]
  sd <- as.numeric(substr(sd,2,nchar(sd)))

  phen <- fread('plink.log',fill=T)
  ncase <- as.numeric(phen$V4[25])
  ncontrol <- as.numeric(phen$V4[24])-ncase
  K <- P <- ncase/(ncase+ncontrol)
  h2.lia <- cvt.h2(K,P,h2)
  h2.sd.lia <- cvt.h2(K,P,sd)
 # print(res)
  print(trait)
  print(h2.lia)
  print(h2.sd.lia)
  #save(res,file='res.lia.RData')
}

for(trait in c( 'HEIGHT','BMI', 'LDL','TG','TC','HDL')){
  setwd('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/')
  dir.create(trait)
  print(trait)
  setwd(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',trait))
  # load('lder.assoc.RData')
  ldpath <- '/home/songs/genomic_control/ukb/20210204/ukb/LD/'
  ldpath.shrink <- '/home/songs/genomic_control/ukb/20210204/ukb/LD.shrink/'
  
  # res <- ld.all(assoc,assoc$N[1],ldpath=ldpath,ldpath.shrink=ldpath.shrink,method='both',
  #               type='jack',rough=F,twostage=T,cores=40)
  res <- fread(paste0('/home/songs/genomic_control/yale_cluster/210428_ukbb_traits/',
                      trait,'/hess/result/step2.log'))
  res <- as.character(res[17,1])
  h2 <- as.numeric(strsplit(res,' ')[[1]][5])
  sd <- strsplit(res,' ')[[1]][6]
  sd <- strsplit(sd,')')[[1]][1]
  sd <- as.numeric(substr(sd,2,nchar(sd)))
  

  # print(res)
  print(trait)
  print(h2)
  print(sd)
  #save(res,file='res.lia.RData')
}


