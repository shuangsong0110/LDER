get.res <- function(x1,lam1,n.gwas,a,rough,twostage){
  m <- length(x1)
  if(twostage){
    newGC <-  1+n.gwas*calH2.new1(x1, lam1,   N=n.gwas,a=NULL, rough=F)$a
    #print(newGC)
    if(newGC<1.05){
      s2thld <- 3.5*sqrt(m/n.gwas)+2.5 
    }else{
      s2thld <- quantile(x1^2, 1-5/m)
    }
    idx <- (x1^2<=s2thld)
    temp1 <-  calH2.new1(x1[idx], lam1[idx],   N=n.gwas,a=a, rough=rough)$a
    a <- max(temp1,0)
    h2 <- calH2.new1(x1, lam1,   N=n.gwas, a=max(temp1,0), rough=rough)$h2
  }else{
    h2 <- calH2.new1(x1, lam1,   N=n.gwas, a=a, rough=rough)$h2
    a <- calH2.new1(x1, lam1,   N=n.gwas, a=a, rough=rough)$a
  }
  return(list(a=a,h2=h2))
}


get.res.ldsc <- function(z1,ldsc1,n.gwas,a,twostage){
  m <- length(z1)
  s1thld <- 30
  if(twostage){
    a <-  calH2(z1[(z1^2)<s1thld], ldsc1[(z1^2)<s1thld],   N=n.gwas,a=NULL)$a
    h2 <- calH2(z1, ldsc1,   N=n.gwas,a=a)$h2
  }else{
    h2 <-  calH2(z1, ldsc1,   N=n.gwas,a=a)$h2
    a <- calH2(z1, ldsc1,   N=n.gwas,a=a)$a
  }
  return(list(a=a,h2=h2))
}
