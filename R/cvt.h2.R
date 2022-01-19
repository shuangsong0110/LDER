cvt.h2 <- function(K,P,h2){
  # K: prevalence
  # P: proportion of cases
  zv <- dnorm(qnorm(K))
  h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
  return(h2_liab)
}
