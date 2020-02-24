ve_fx <- function(data, freq) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  #number of categories
  num_cat <- sapply(data, function(x) length(unique(x)))
  
  #frequency tables
  rel.freq <- freq/sum(freq[,1])
  ln.freq <- log(rel.freq)
  ln.freq[ln.freq == -Inf] <- 0
  
  #entropy
  entropy_matrix <- rel.freq * ln.freq
  entropy<- - colSums(entropy_matrix)
  norm_entropy <- entropy/log(num_cat)
  norm_entropy <- ifelse(is.nan(norm_entropy),0,norm_entropy)
  
  agreement <- vector(mode="numeric", length=s)
  ve <- matrix(data=0,nrow=r,ncol=r)
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        if (data[i,k] == data[j,k]) {
          agreement[k] <- norm_entropy[k]
        }
        else {
          agreement[k] <- 0
        }
      }
      if (i == j) {
        ve[i,j] <- 0
      }
      else {
        ve[i,j] <- 1-1/s*(sum(agreement))
        ve[j,i] <- ve[i,j]
      }
    }
  }
  return(ve)
}