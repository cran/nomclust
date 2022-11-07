GetFreqRel <- function(freq_table) {

  freq_rel <- freq_table/sum(freq_table[,1])
  
  return(freq_rel)
}


GetFreqRel2 <- function(freq_table) {
  
  freq_rel2 <-  freq_table * (freq_table -1) / r / (r-1)
  
  return(freq_rel2)
}
