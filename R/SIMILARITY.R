SIMILARITY <- function(data, measure, freq.table, wt) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  num_cat <- sapply(data, function(x) length(unique(x)))
  num_cat_sum <- sum(num_cat)
  freq.rel <- freq.table/sum(freq.table[,1])
  freq.rel.simple2 <- freq.rel^2
  freq.rel2 <- freq.table * (freq.table -1) / r / (r-1)
  #wt <- rep(1, s) # IF variable weights are turned-off
  
  if (measure == "eskin") {
    prox_matrix <- matrix(eskin_cpp(r, s, num_cat, as.vector(t(t(data))), wt, sum(wt)), ncol=r, nrow = r)
  } else if (measure == "anderberg") {
    prox_matrix <- matrix(anderberg_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), nrow(freq.rel), num_cat, wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "burnaby") {
    prox_matrix <- matrix(burnaby_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "gambaryan") {
    prox_matrix <- matrix(gambaryan_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), nrow(freq.rel), num_cat_sum, wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "goodall1") {
    prox_matrix <- matrix(good1_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), as.vector(t(t(freq.rel2))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "goodall2") {
    prox_matrix <- matrix(good2_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), as.vector(t(t(freq.rel2))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "goodall3") {
    prox_matrix <- matrix(good3_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel2))), nrow(freq.rel2), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "goodall4") {
    prox_matrix <- matrix(good4_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel2))), nrow(freq.rel2), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "iof") {
    prox_matrix <- matrix(iof_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.table))), nrow(freq.table), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "lin") {
    prox_matrix <- matrix(lin_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
    prox_matrix[prox_matrix == -Inf] <- max(prox_matrix) + 1
  } else if (measure == "lin1") {
    prox_matrix <- matrix(lin1_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.rel))), nrow(freq.rel), wt, sum(wt)), ncol=r, nrow=r)
    prox_matrix[prox_matrix == -Inf] <- max(prox_matrix) + 1
  } else if (measure == "of") {
    prox_matrix <- matrix(of_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.table))), nrow(freq.table), sum(freq.table[,1]), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "sm") {
    prox_matrix <- matrix(sm_cpp(r, s, as.vector(t(t(data))), wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "smirnov") {
    prox_matrix <- matrix(smirnov_cpp(r, s, as.vector(t(t(data))), as.vector(t(t(freq.table))), nrow(freq.table), num_cat_sum, wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "ve") {
    #number of categories
    ln.freq <- log(freq.rel)
    ln.freq[ln.freq == -Inf] <- 0
    #entropy
    entropy_matrix <- freq.rel * ln.freq
    entropy<- - colSums(entropy_matrix)
    norm_entropy <- entropy/log(num_cat)
    norm_entropy <- ifelse(is.nan(norm_entropy),0,norm_entropy)
    prox_matrix <- matrix(ve_cpp(r, s, as.vector(t(t(data))), norm_entropy, wt, sum(wt)), ncol=r, nrow=r)
  } else if (measure == "vm") {
    #gini coefficient
    sum_rel2.freq <- colSums(freq.rel.simple2)
    gini <- 1- sum_rel2.freq
    norm_gini <- gini*num_cat/(num_cat-1)
    norm_gini <- ifelse(is.nan(norm_gini),0,norm_gini)
    prox_matrix <- matrix(vm_cpp(r, s, as.vector(t(t(data))), norm_gini, wt, sum(wt)), ncol=r, nrow=r)
  } else{
    stop("Invalid name of the similarity measure.")
  }
    
  return(prox_matrix)
  
}