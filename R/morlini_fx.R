
morlini_fx <- function(data, freq = NULL) {
  
  s <- ncol(data)
  num_cat <- sapply(data, function(x) length(unique(x)))
  
  rnames <- row.names(data)
  
  #with dummies
  #data_dummy <- dummy.data.frame(data, dummy.classes ="ALL",omit.constants = F)
  
  # dummy transformation of a data.frame
  data_dummy <- data.frame(row.names = row.names(data))
  num_cat <- sapply(data, function(x) length(unique(x)))
  for (i in 1:length(num_cat)) {
    variable_set <- data.frame()
    for (j in 1:num_cat[i]) {
      for (k in (1:nrow(data))) {
        variable_set[k,j] <- ifelse(data[k,i] == j, 1, 0)
      }
    }
    data_dummy <- cbind(data_dummy, variable_set)
  }
  
  n <- nrow(data_dummy)
  hs <- ncol(data_dummy)
  
  nsv <- sapply(data_dummy, sum)
  fsv2 <- log(1/(nsv/n)^2)
  
  E <- matrix(data=0,nrow=n,ncol=n)
  agreement <- vector(mode="numeric", length=hs)
  
  #computation of Eij
  for (i in 1:(n-1)) {
    for (j in (1+i):n) {
      for (k in 1:hs) {
        if (data_dummy[i,k] == 1 & data_dummy[j,k] == 1) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- 0
        }
      }
      E[i,j] <- fsv2 %*% agreement
      E[j,i] <- E[i,j]
    }
  }

#computation of Fij
  cum <- cumsum(num_cat)
  F <- matrix(data=0,nrow=n,ncol=n)
  row.names(F) <- rnames
  
  for (i in 1:(n-1)) {
    for (j in (1+i):n) {
      v <- 0
      agreement <- vector(mode="numeric", length=hs)
      for (k in 1:s) {
        for (t in (v+1):cum[k]) {
          if (data_dummy[i,t] == 0 & data_dummy[j,t] == 1) {
            agreement[(v+1):cum[k]] <- 1
          }
        }
        v <- cum[k]      
      }
      F[i,j] <- t(agreement) %*% fsv2
      F[j,i] <- F[i,j]
    }
  }
   morlini <- 1 - E/(E+F)
   morlini <- ifelse(is.nan(morlini),0,morlini)
  return(morlini)
}
  


  
  