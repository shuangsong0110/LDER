bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm=T)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}
