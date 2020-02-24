lin_fx <- function(data, freq) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  freq.rel <- freq/sum(freq[,1])
  
  agreement <- vector(mode="numeric", length=s)
  lin <- matrix(data=0,nrow=r,ncol=r)
  weights <- vector(mode="numeric", length=s)
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        d <- data[j,k]
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 2*log(freq.rel[c,k])
        }
        else {
          agreement[k] <- 2*log(freq.rel[c,k] + freq.rel[d,k])
        }
        weights[k] <- log(freq.rel[c,k]) + log(freq.rel[d,k])
      }
      if (i == j) {
        lin[i,j] <- 0
      }
      else {
        lin[i,j] <- 1/(1/sum(weights)*(sum(agreement))) - 1
        lin[j,i] <- lin[i,j]
      }
    }
  }
  lin[lin == -Inf] <- max(lin) + 1
  return(lin)
}