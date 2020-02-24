sm_fx <- function(data, freq = NULL) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  agreement <- vector(mode="numeric", length=s)
  sm <- matrix(data=0,nrow=r,ncol=r)
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- 0
        }
      }
      sm[i,j] <- 1-1/s*(sum(agreement))
      sm[j,i] <- sm[i,j]
    }
  }
  return(sm)
}
