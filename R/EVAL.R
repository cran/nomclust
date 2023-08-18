EVAL <- function(M, clusters, diss){

    # determination of minimal and maximal number of clusters
  clu_high <- dim(M[[length(M)]])[2]
  clu_low <- clu_high - length(M) + 2
  num_clu <-  clu_high - clu_low + 1
  
  #number of variables
  m <- dim(M[[clu_low]])[3]

  #number of elements
  n <- sum(M[[clu_low]][,,1])
  
  # number of categories
  num_cat <- vector(mode="numeric", length=m)
  for (i in 1:m) {
    num_cat[i] <- sum(M[[1]][,,i]!=0)
  }
  
  names <- vector(mode="character", length=num_clu)
  BIC <- rep(NA, num_clu)
  AIC <- rep(NA, num_clu)
  WCE <- rep(NA, num_clu)
  WCM <- rep(NA, num_clu)
  nWCE <- rep(NA, num_clu)
  nWCM <- rep(NA, num_clu)
  PSFE <- rep(NA, num_clu)
  PSFM <- rep(NA, num_clu)
  increment <- rep(NA, num_clu)
  BK <- rep(NA, num_clu)
  SI <- rep(NA, num_clu)
  CU <- rep(NA, num_clu)
  CU_1 <- rep(NA, num_clu)
  CI <- rep(NA, num_clu)
  HE <- rep(NA, num_clu)
  HM <- rep(NA, num_clu)
  nHE <- rep(NA, num_clu)
  nHM <- rep(NA, num_clu)
  DI <- rep(NA, num_clu)
  
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
    nWCE[k] = sum(H_g)/n # non-normalized
    nWCM[k] = sum(G_g)/n # non-normalized
    PSFE[k] <- ((n-k) * (nWCE[1] - nWCE[k]))/((k-1)*nWCE[k]) # BASED ON NON-NORMALIZED WCE
    PSFE[1] <- NA
    PSFM[k] <- ((n-k) * (nWCM[1] - nWCM[k]))/((k-1)*nWCM[k]) # BASED ON NON-NORMALIZED WCM
    PSFM[1] <- NA
    increment[k] <- sum(H_g_norm)
    CU[k] <- (nWCM[1] - nWCM[k]) / k
    CU[1] <- NA
    CI[k] <- (nWCE[1] - nWCE[k]) / k
    CI[1] <- NA
    }
  
  
  # Hartigan indices
  for (k in 1:clu_high) {
    #HE[k] <- (WCE[k] / WCE[k+1] - 1) * (n - k - 1)
    #HM[k] <- (WCM[k] / WCM[k+1] - 1) * (n - k - 1)
    HE[k] <- (nWCE[k] / nWCE[k+1] - 1) * (n - k - 1)
    HM[k] <- (nWCM[k] / nWCM[k+1] - 1) * (n - k - 1)
  }
  

  # BK index
  for (k in clu_low:clu_high) {
    BK[k] <- (increment[k-1]-increment[k]) - (increment[k]-increment[k+1])
    BK[1] <- NA
  }
  
  # CU
 # CU <- rep(NA, clu_high)
  # for (k in 2:clu_high){
  #   # cluster sizes
  #   clu_size <- vector(mode="numeric", length = k)
  #   for (i in 1:k) {
  #     clu_size[i] <- sum(M[[k]][,i,1])
  #   }
  #   # relative sizes
  #   p_clu <- clu_size / n 
  #   
  #   # unconditional probabilities
  #   p_uncond <- sum((M[[1]][,1,] / sum(clu_size))^2)
  #   
  #   # conditional probabilities
  #   p_cond <- vector(mode="numeric", length = k)
  #   nominator <- 0
  #   for (i in 1:k) {
  #     p_cond[i] <- sum((M[[k]][,i,] / clu_size[i])^2)
  #     output <- p_clu[i] * (p_cond[i] - p_uncond)
  #     nominator <- nominator + output
  #   }
  #   
  #   CU[k] <- nominator / k
  # }
  
  
  # dissimilarity matrix is in the input
  if (is.null(diss) == FALSE) {
    for (k in clu_low:clu_high) {
      SI[k] <- summary(silhouette(clusters[,k-1], diss))$avg.width
    }
    for (k in clu_low:clu_high) {
      DI[k] <- dunn(distance = diss, clusters[,k-1], Data = NULL)
    }
    
  }
  
  if (is.null(diss) == FALSE) {
  
  # list of coefficients in an output
  coef <- list(names = names, WCM = WCM, WCE = WCE, PSFM = PSFM, PSFE = PSFE, BIC = BIC, AIC = AIC, BK = BK, SI = SI, DI = DI,
               CU = CU, CI = CI, HE = HE, HM = HM)

  clusters <- c(1,clu_low:clu_high) # a list of examined cluster solutions
  PSFM_opt <- as.integer(clusters[which.max(coef$PSFM)])
  PSFE_opt <- as.integer(clusters[which.max(coef$PSFE)])
  BK_opt <- as.integer(clusters[which.max(coef$BK)])
  BIC_opt <- as.integer(clusters[which.min(coef$BIC)])
  AIC_opt <- as.integer(clusters[which.min(coef$AIC)])
  SI_opt <- as.integer(clusters[which.max(coef$SI)])
  CU_opt <- as.integer(clusters[which.max(coef$CU)])
  CI_opt <- as.integer(clusters[which.max(coef$CI)])
  HE_opt <- as.integer(clusters[which.min(coef$HE)])
  HM_opt <- as.integer(clusters[which.min(coef$HM)])
  DI_opt <- as.integer(clusters[which.max(coef$DI)])

  optimal <- list(PSFM = PSFM_opt, PSFE = PSFE_opt, BIC = BIC_opt, AIC = AIC_opt, BK = BK_opt, SI = SI_opt, DI = DI_opt, CU = CU_opt, CI = CI_opt, HE = HE_opt, HM = HM_opt)
  } else {
    # list of coefficients in an output
    coef <- list(names = names, WCM = WCM, WCE = WCE, PSFM = PSFM, PSFE = PSFE, BIC = BIC, AIC = AIC, BK = BK, CU = CU, CI = CI, HE = HE, HM = HM)
    
    clusters <- c(1,clu_low:clu_high) # a list of examined cluster solutions
    PSFM_opt <- as.integer(clusters[which.max(coef$PSFM)])
    PSFE_opt <- as.integer(clusters[which.max(coef$PSFE)])
    BK_opt <- as.integer(clusters[which.max(coef$BK)])
    BIC_opt <- as.integer(clusters[which.min(coef$BIC)])
    AIC_opt <- as.integer(clusters[which.min(coef$AIC)])
    CU_opt <- as.integer(clusters[which.max(coef$CU)])
    CI_opt <- as.integer(clusters[which.max(coef$CI)])
    HE_opt <- as.integer(clusters[which.min(coef$HE)])
    HM_opt <- as.integer(clusters[which.min(coef$HM)])
    
    optimal <- list(PSFM = PSFM_opt, PSFE = PSFE_opt, BIC = BIC_opt, AIC = AIC_opt, BK = BK_opt, CU = CU_opt, CI = CI_opt, HE = HE_opt, HM = HM_opt)
  }
  
  output <- list(coef, optimal) 
  return(output)
}


