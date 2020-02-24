vm_fx <- function(data, freq) {

  r <- nrow(data)
  s <- ncol(data)
  
  #number of categories
  num_cat <- sapply(data, function(x) length(unique(x)))

  #frequency tables
  rel.freq <- freq/sum(freq[,1])
  rel2.freq <- rel.freq^2

  #gini coefficient
  sum_rel2.freq <- colSums(rel2.freq)
  gini <- 1- sum_rel2.freq
  norm_gini <- gini*num_cat/(num_cat-1)
  norm_gini <- ifelse(is.nan(norm_gini),0,norm_gini)
  
  agreement <- vector(mode="numeric", length=s)
  vm <- matrix(data=0,nrow=r,ncol=r)

  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        if (data[i,k] == data[j,k]) {
          agreement[k] <- norm_gini[k]
       }
       else {
          agreement[k] <- 0
        }
      }
      if (i == j) {
        vm[i,j] <- 0
      }
      else {
        vm[i,j] <- 1-1/s*(sum(agreement))
        vm[j,i] <- vm[i,j]
      }
    }
  }
  return(vm)
}