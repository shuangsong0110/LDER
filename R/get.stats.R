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
