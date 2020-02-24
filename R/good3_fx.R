good3_fx <- function(data, freq) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  freq.rel <- freq/sum(freq[,1])
  
  agreement <- vector(mode="numeric", length=s)
  good3 <- matrix(data=0,nrow=r,ncol=r)
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 1 - freq.rel[c,k]^2
        }
        else {
          agreement[k] <- 0
        }
      }
      if (i == j) {
        good3[i,j] <- 0
      }
      else {
        good3[i,j] <- 1-1/s*(sum(agreement))
        good3[j,i] <- good3[i,j]
      }
    }
  }
  return(good3)
}