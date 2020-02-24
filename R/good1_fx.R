good1_fx <- function(data, freq) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  freq.rel <- freq/sum(freq[,1])
  freq.rel2 <- freq.rel^2
  
  agreement <- vector(mode="numeric", length=s)
  good1 <- matrix(data=0,nrow=r,ncol=r)
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        if (data[i,k] == data[j,k]) {
          logic <- t(freq.rel[,k] <= freq.rel[c,k])
          agreement[k] <- 1 - sum(freq.rel2[,k] * logic)
        }
        else {
          agreement[k] <- 0
        }
      }
      if (i == j) {
        good1[i,j] <- 0
      }
      else {
        good1[i,j] <- 1-1/s*(sum(agreement))
        good1[j,i] <- good1[i,j]
      }
    }
  }
  return(good1)
}