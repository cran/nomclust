SIMILARITY <- function(data, measure, freq.table) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  wt <- rep(1, s) # IF variable weigths are turned-off
  
  if (measure == "eskin") {
    num_cat <- sapply(data, function(x) length(unique(x)))
    prox_matrix <- matrix(eskin_cpp(r, s, num_cat, as.vector(t(t(data))), wt, sum(wt)), ncol=r, nrow = r)
  } else if (measure == "good1") {
    freq.rel <- freq.table/sum(freq.table[,1])
    freq.rel2 <- freq.rel^2
    prox_matrix <- matrix(good1_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), as.vector(t(t(freq.rel2))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "good2") {
    freq.rel <- freq.table/sum(freq.table[,1])
    freq.rel2 <- freq.rel^2
    prox_matrix <- matrix(good2_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), as.vector(t(t(freq.rel2))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "good3") {
    freq.rel <- freq.table/sum(freq.table[,1])
    prox_matrix <- matrix(good3_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "good4") {
    freq.rel <- freq.table/sum(freq.table[,1])
    prox_matrix <- matrix(good4_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "iof") {
    freq <- freq.table
    prox_matrix <- matrix(iof_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq))), nrow(freq), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "lin") {
    freq.rel <- freq.table/sum(freq.table[,1])
    prox_matrix <- matrix(lin_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
    prox_matrix[prox_matrix == -Inf] <- max(prox_matrix) + 1
  } else if (measure == "lin1") {
    freq.rel <- freq.table/sum(freq.table[,1])
    prox_matrix <- matrix(lin1_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
    prox_matrix[prox_matrix == -Inf] <- max(prox_matrix) + 1
  } else if (measure == "of") {
    n <- sum(freq.table[,1])
    prox_matrix <- matrix(of_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.table))), nrow(freq.table), n, wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "sm") {
    prox_matrix <- matrix(sm_cpp(r, s, as.vector(t(t(data))), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "ve") {
    #number of categories
    num_cat <- sapply(data, function(x) length(unique(x)))
    #frequency tables
    rel.freq <- freq.table/sum(freq.table[,1])
    ln.freq <- log(rel.freq)
    ln.freq[ln.freq == -Inf] <- 0
    #entropy
    entropy_matrix <- rel.freq * ln.freq
    entropy<- - colSums(entropy_matrix)
    norm_entropy <- entropy/log(num_cat)
    norm_entropy <- ifelse(is.nan(norm_entropy),0,norm_entropy)
    prox_matrix <- matrix(ve_cpp(r, s, as.vector(t(t(data))), norm_entropy, wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "vm") {
    #number of categories
    num_cat <- sapply(data, function(x) length(unique(x)))
    #frequency tables
    rel.freq <- freq.table/sum(freq.table[,1])
    rel2.freq <- rel.freq^2
    #gini coefficient
    sum_rel2.freq <- colSums(rel2.freq)
    gini <- 1- sum_rel2.freq
    norm_gini <- gini*num_cat/(num_cat-1)
    norm_gini <- ifelse(is.nan(norm_gini),0,norm_gini)
    prox_matrix <- matrix(vm_cpp(r, s, as.vector(t(t(data))), norm_gini, wt, sum(wt)), ncol=r, nrow=r)
  } else{
    stop("Invalid name of the similarity measure.")
  }
    
  return(prox_matrix)
  
}