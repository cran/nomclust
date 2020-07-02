lin1_fx <- function(data, freq) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  freq.rel <- freq/sum(freq[,1])
  #freq.ln <- log(freq.rel)
  #freq.ln[freq.ln == -Inf] <- 0
  
  
  agreement <- vector(mode="numeric", length=s)
  lin1 <- matrix(data=0,nrow=r,ncol=r)
  weights <- vector(mode="numeric", length=s)

  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        d <- data[j,k]
        if (data[i,k] == data[j,k]) {
          logic <- freq.rel[,k] == freq.rel[c,k]
          agreement[k] <- log(sum(logic * freq.rel[,k]))
          weights[k] <- log(freq.rel[c,k]) + log(freq.rel[d,k])
        }
        else {
          if (freq.rel[c,k] >= freq.rel[d,k]) {
            logic <- freq.rel[,k] >= freq.rel[d,k] & freq.rel[,k] <= freq.rel[c,k]
            agreement[k] <- 2*log(sum(logic * freq.rel[,k]))
            weights[k] <- log(freq.rel[c,k]) + log(freq.rel[d,k])
          }
          else {
            logic <- freq.rel[,k] >= freq.rel[c,k] & freq.rel[,k] <= freq.rel[d,k]
            agreement[k] <- 2*log(sum(logic * freq.rel[,k]))
            weights[k] <- log(freq.rel[c,k]) + log(freq.rel[d,k])
          }
        }
      }
      lin1[i,j] <- 1/(1/sum(weights)*(sum(agreement))) - 1
      lin1[j,i] <- lin1[i,j]
    }
  }
  lin1[lin1 == -Inf] <- max(lin1) + 1
  return(lin1)
}