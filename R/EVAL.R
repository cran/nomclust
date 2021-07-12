EVAL <- function(M, clusters, diss){

    # determination of minimal and maximal number of clusters
  clu_high <- dim(M[[length(M)]])[2]
  clu_low <- clu_high - length(M) + 2
  num_clu <-  clu_high - clu_low + 1
  
  #number of variables
  m <- dim(M[[clu_low]])[3]

  #number of elements
  n <- sum(M[[clu_low]][,,1])
  
  num_cat <- vector(mode="numeric", length=m)
  for (i in 1:m) {
    num_cat[i] <- sum(M[[1]][,,i]!=0)
  }
  
  names <- vector(mode="character", length=num_clu)
  BIC <- vector(mode="numeric", length=num_clu)
  AIC <- vector(mode="numeric", length=num_clu)
  WCE <- vector(mode="numeric", length=num_clu)
  WCM <- vector(mode="numeric", length=num_clu)
  nWCE <- vector(mode="numeric", length=num_clu)
  nWCM <- vector(mode="numeric", length=num_clu)
  PSFE <- vector(mode="numeric", length=num_clu)
  PSFM <- vector(mode="numeric", length=num_clu)
  increment <- vector(mode="numeric", length=num_clu)
  BK <- vector(mode="numeric", length=num_clu)
  SI <- vector(mode="numeric", length=num_clu)
  
    for (k in c(1,clu_low:clu_high)) {
    matrix <- M[[k]]
    
    H_g <- vector(mode="numeric", length=k)
    H_g_norm <- vector(mode="numeric", length=k)
    G_g <- vector(mode="numeric", length=k)
    G_g_norm <- vector(mode="numeric", length=k)
    
    penalty <- sum(num_cat - 1)

    for (g in 1:k) {
      H_gc <- vector(mode="numeric", length=m)
      G_gc <- vector(mode="numeric", length=m)
      H_gc_norm <- vector(mode="numeric", length=m)
      G_gc_norm <- vector(mode="numeric", length=m)
      n_g <- sum(matrix[,g,1])
      
      for (c in 1:m) {
        step <- 0
        step <- matrix[,g,c]
        
        K <- num_cat[c]

        temp_H <- vector(mode="numeric", length=K)
        temp_G <- vector(mode="numeric", length=K)
        for (u in 1:K) {
          if (step[u] == 0) {
            temp_H[u] <- 0
            temp_G[u] <- 0
          }
          else {
            temp_H[u] <- ((step[u]/n_g)%*%log(step[u]/n_g))
            temp_G[u] <- ((step[u]/n_g)^2)
          }
        }
        
        H_gc[c] <- -sum(temp_H)
        H_gc[c] <- ifelse(is.nan(H_gc[c]),0,H_gc[c])
        G_gc[c] <- 1-sum(temp_G)
        H_gc_norm[c] <- H_gc[c]/log(K)
        G_gc_norm[c] <- G_gc[c]*K/(K-1)
      }
      
      H_g[g] <- sum(H_gc) * n_g # for BIC, AIC
      G_g[g] <- sum(G_gc) * n_g # for BIC, AIC
      H_g_norm[g] <- sum(H_gc_norm) * n_g /n  # for BK, WCE
      G_g_norm[g] <- sum(G_gc_norm) * n_g /n # for WCM
      
    }
    names[k] <- paste("clu",k, sep = "_")
    BIC[k] <- -2 * -(sum(H_g)) + k * penalty * log(n)
    AIC[k] = -2 * -(sum(H_g)) + 2*k * penalty
    WCE[k] = sum(H_g_norm)/m # normalized
    WCM[k] = sum(G_g_norm)/m # normalized
    nWCE[k] = sum(H_g)/n/m # non-normalized
    nWCM[k] = sum(G_g)/n/m # non-normalized
    PSFE[k] <- ((n-k) * (nWCE[1] - nWCE[k]))/((k-1)*nWCE[k]) # BASED ON NON-NORMALIZED WCE
    PSFE[1] <- NA
    PSFM[k] <- ((n-k) * (nWCM[1] - nWCM[k]))/((k-1)*nWCM[k]) # BASED ON NON-NORMALIZED WCM
    PSFM[1] <- NA
    increment[k] <- sum(H_g_norm)
  }

  for (k in clu_low:clu_high) {
    BK[k] <- (increment[k-1]-increment[k]) - (increment[k]-increment[k+1])
    BK[1] <- NA
  }
  
  # dissimilarity matrix is in the input
  if (is.null(diss) == FALSE) {
    for (k in clu_low:clu_high) {
      SI[k] <- summary(silhouette(clusters[,k-1], diss))$avg.width
      SI[1] <- NA
    }
  }
  
  if (is.null(diss) == FALSE) {
  
  # list of coefficients in an output
  coef <- list(names = names, WCM = WCM, WCE = WCE, PSFM = PSFM, PSFE = PSFE, BIC = BIC, AIC = AIC, BK = BK, SI = SI)

  clusters <- c(1,clu_low:clu_high) # a list of examined cluster solutions
  PSFM_opt <- as.integer(clusters[which.max(coef$PSFM)])
  PSFE_opt <- as.integer(clusters[which.max(coef$PSFE)])
  BK_opt <- as.integer(clusters[which.max(coef$BK)])
  BIC_opt <- as.integer(clusters[which.min(coef$BIC)])
  AIC_opt <- as.integer(clusters[which.min(coef$AIC)])
  SI_opt <- as.integer(clusters[which.max(coef$SI)])

  optimal <- list(PSFM = PSFM_opt, PSFE = PSFE_opt, BIC = BIC_opt, AIC = AIC_opt, BK = BK_opt, SI = SI_opt)
  } else {
    # list of coefficients in an output
    coef <- list(names = names, WCM = WCM, WCE = WCE, PSFM = PSFM, PSFE = PSFE, BIC = BIC, AIC = AIC, BK = BK)
    
    clusters <- c(1,clu_low:clu_high) # a list of examined cluster solutions
    PSFM_opt <- as.integer(clusters[which.max(coef$PSFM)])
    PSFE_opt <- as.integer(clusters[which.max(coef$PSFE)])
    BK_opt <- as.integer(clusters[which.max(coef$BK)])
    BIC_opt <- as.integer(clusters[which.min(coef$BIC)])
    AIC_opt <- as.integer(clusters[which.min(coef$AIC)])
    
    optimal <- list(PSFM = PSFM_opt, PSFE = PSFE_opt, BIC = BIC_opt, AIC = AIC_opt, BK = BK_opt)
  }
  
  output <- list(coef, optimal) 
  return(output)
}


