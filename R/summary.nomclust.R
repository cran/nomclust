summary.nomclust <- function(object, ...){

  if (names(object)[2] == "eval"){ # nomclust object and nomprox with data
  
    cat("\nSizes of the created clusters:\n")
    for (i in 1:length(object$mem)) {
      cat("\n", i+1," clusters:", sep = "")
      print(table(object[[1]][i]))
    }
    
    cat("\nOptimal number of clusters based on the evaluation criteria:\n")
    print(as.data.frame(object$opt))
    
    cat("\nAgglomerative coefficient:", object$dend$ac, "\n")
  }
  
  if (names(object)[2] == "opt"){ # evalclust
    cat("\nOptimal number of clusters based on the evaluation criteria:\n")
    print(as.data.frame(object$opt))
  }
  
  
  if (names(object)[2] == "dend"){ # nomprox without data
    cat("\nSizes of the created clusters:\n")
    for (i in 1:length(object$mem)) {
      cat("\n", i+1," clusters:", sep = "")
      print(table(object[[1]][i]))
    }
    cat("\nAgglomerative coefficient:", object$dend$ac, "\n")
  }
}
