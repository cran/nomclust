#' Variable Mutability (VM) measure
#'
#' @description A function for calculation of a proximity (dissimilarity) matrix based on the VM similarity measure.
#' \cr                                                           
#'  
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#' 
#' @return The function returns a dissimilarity matrix of the size \code{n x n}, where \code{n} is the number of objects in the original dataset in the argument \code{data}.
#' \cr
#' 
#' @details The Variable Mutability similarity measure was introduced in (Sulc and Rezankova, 2019).
#' It treats the similarity between two categories based on the within-cluster variability expressed by the normalized mutability. The measure assigns higher weights to rarer categories. 
#' 
#' @references
#' Sulc Z. and Rezankova H. (2019). Comparison of Similarity Measures for Categorical Data in Hierarchical Clustering. Journal of Classification. 2019, 35(1), p. 58-72. DOI: 10.1007/s00357-019-09317-5.
#'
#' @seealso
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{good1}},
#' \code{\link[nomclust]{good2}},
#' \code{\link[nomclust]{good3}},
#' \code{\link[nomclust]{good4}},
#' \code{\link[nomclust]{iof}},
#' \code{\link[nomclust]{lin}},
#' \code{\link[nomclust]{lin1}},
#' \code{\link[nomclust]{morlini}}, 
#' \code{\link[nomclust]{of}},
#' \code{\link[nomclust]{sm}},
#' \code{\link[nomclust]{ve}}.
#'
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' #sample data
#' data(data20)
#' 
#' # dissimilarity matrix calculation
#' prox.vm <- vm(data20)
#' 
#' @export 


vm <- function(data) {

  r <- nrow(data)
  s <- ncol(data)
  
  rnames <- row.names(data)
  
  # recoding everything to factors and then to numeric values
  indx <- sapply(data, is.factor)
  data[!indx] <- sapply(data[!indx], function(x) as.factor(x))
  data <- as.data.frame(unclass(data))
  data <- sapply(data, function(x) as.numeric(x))
  data <- as.data.frame(data)


  #number of categories
  num_cat <- sapply(data, function(x) length(unique(x)))

  #frequency tables
  abs.freq <- freq.abs(data)
  rel.freq <- abs.freq/r
  rel2.freq <- rel.freq^2

  #gini coefficient
  sum_rel2.freq <- colSums(rel2.freq)
  gini <- 1- sum_rel2.freq
  norm_gini <- gini*num_cat/(num_cat-1)
  norm_gini <- ifelse(is.nan(norm_gini),0,norm_gini)
  
  agreement <- vector(mode="numeric", length=s)
  vm <- matrix(data=0,nrow=r,ncol=r)
  row.names(vm) <- rnames
  
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