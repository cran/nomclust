#' Variable Mutability measure
#' 
#' @description The Variable Mutability similarity measure was introduced in (Sulc and Rezankova, 2015).
#' It treats similarity between two categories according to within-cluster variability expressed by the Gini coefficient (mutability).
#' The novel similarity measures praise more the match of two categories in a variable with high variability, because it is rarer,
#' than the match in a low-variability variable.
#' Hierarchical clustering methods require a proximity (dissimilarity) matrix instead of a similarity matrix as
#' an entry for the analysis; therefore, dissimilarity \code{D} is computed from similarity \code{S} according the equation
#' \code{1/S-1}.\cr
#' \cr                                                           
#'  
#' @param data data frame or matrix with cases in rows and variables in colums. Cases are characterized by nominal (categorical) variables coded as numbers.
#' 
#' @return Function returns a matrix of the size \code{n x n}, where \code{n} is the number of objects in original data. The matrix contains proximities
#' between all pairs of objects. It can be used in hierarchical cluster analyses (HCA), e.g. in \code{\link[cluster]{agnes}}.
#' \cr
#' @references
#' Sulc, Z. and Rezankova H. (2015). Novel similarity measures for categorical data based on mutability and entropy.
#'  Conference of the International Federation of Classification Societies. Bologna: Ospitalia, p. 209.
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
#' # Creation of proximity matrix
#' prox_vm <- vm(data20)
#' 
#' @export 


vm <- function(data) {

  r <- nrow(data)
  s <- ncol(data)
  
  #recoding variables
  num_var <- ncol(data)
  num_row <- nrow(data)
  data2 <- matrix(data = 0, nrow = num_row, ncol = num_var)
  for (k in 1:num_var) {
    categories <- unique(data[, k])
    cat_new <- 1:length(categories)
    for (l in 1:length(categories)) {
      for (i in 1:num_row) {
        if (data[i, k] == categories[l]) {
          data2[i, k] <- cat_new[l]
        }
      }
    }
  }
  data <- data.frame(data2)


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