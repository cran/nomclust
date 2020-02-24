eskin_fx <- function(data, freq = NULL) {
  
  r <- nrow(data)
  s <- ncol(data)

  num_cat <- sapply(data, function(x) length(unique(x)))
  
  agreement <- vector(mode="numeric", length=s)
  eskin <- matrix(data=0,nrow=r,ncol=r)
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- num_cat[k]^2/(num_cat[k]^2 + 2)
        }
      }
      eskin[i,j] <- 1/(1/s*(sum(agreement))) - 1
      eskin[j,i] <- eskin[i,j]
    }
  }
  return(eskin)
}
