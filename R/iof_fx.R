iof_fx <- function(data, freq) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  agreement <- vector(mode="numeric", length=s)
  iof <- matrix(data=0,nrow=r,ncol=r)
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        d <- data[j,k]
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- 1/(1+log(freq[c,k])*log(freq[d,k]))
        }
      }
      iof[i,j] <- 1/(1/s*(sum(agreement)))-1
      iof[j,i] <- iof[i,j]
    }
  }
  return(iof)
}



