of_fx <- function(data, freq) {
  
  r <- nrow(data)
  s <- ncol(data)
  n <- sum(freq[,1])

  agreement <- vector(mode="numeric", length=s)
  of <- matrix(data=0,nrow=r,ncol=r)
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        d <- data[j,k]
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- 1/(1+log(n/freq[c,k])*log(n/freq[d,k]))
        }
      }
      of[i,j] <- 1/(1/s*(sum(agreement)))-1
      of[j,i] <- of[i,j]
    }
  }
  return(of)
}

