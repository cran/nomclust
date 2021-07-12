print.nomclust <- function(x, ...){
  #cat("\n Call: ")
  #print(match.call())
  if (names(x)[2] == "eval"){ # nomclust object and nomprox with data
    cat("\nValues of the evaluation criteria:\n")
    print(as.data.frame(x$eval))
    
    cat("\nOptimal number of clusters based on the evaluation criteria:\n")
    print(as.data.frame(x$opt))
    
    cat("\nAgglomerative coefficient:", x$dend$ac, "\n")
  }
  
  if (names(x)[2] == "opt"){ # evalclust
    cat("\nValues of the evaluation criteria:\n")
    print(as.data.frame(x$eval))
    
    cat("\nOptimal number of clusters based on the evaluation criteria:\n")
    print(as.data.frame(x$opt))
  }
  
  
  if (names(x)[2] == "dend"){ # nomprox without data
    cat("\nAgglomerative coefficient:", x$dend$ac, "\n")
  }
   
}
