#' Variable Entropy (VE) Measure
#' 
#' @description A function for calculation of a proximity (dissimilarity) matrix based on the VE similarity measure.
#'  
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#' 
#' @return The function returns an object of class "dist".
#' \cr
#' 
#' @details The Variable Entropy similarity measure was introduced in (Sulc and Rezankova, 2019). It treats
#' the similarity between two categories based on the within-cluster variability expressed by the normalized entropy.
#' The measure assigns higher weights to rare categories.
#' 
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'  \cr
#'  \cr
#' Sulc Z. and Rezankova H. (2019). Comparison of Similarity Measures for Categorical Data in Hierarchical Clustering. Journal of Classification. 2019, 35(1), p. 58-72. DOI: 10.1007/s00357-019-09317-5.
#'  
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
#' \code{\link[nomclust]{of}},
#' \code{\link[nomclust]{sm}},
#' \code{\link[nomclust]{vm}}.
#'
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' # sample data
#' data(data20)
#' 
#' # dissimilarity matrix calculation
#' prox.ve <- ve(data20)
#' 
#' @export 

ve <- function(data) {
  
  # dealing with the missing data
  if (sum(is.na(data)) > 0) {
    stop("The dissimilarity matrix CANNOT be calculated if the 'data' argument contains NA values.")
  }
  
  r <- nrow(data)
  s <- ncol(data)
  
  rnames <- row.names(data)
  
  # recoding everything to factors and then to numeric values
  indx <- sapply(data, is.factor)
  data[!indx] <- lapply(data[!indx], function(x) as.factor(x))
  data <- as.data.frame(sapply(data, function(x) as.numeric(x)))
  
  
  #number of categories
  num_cat <- sapply(data, function(x) length(unique(x)))
  
  #frequency tables
  abs.freq <- freq.abs(data)
  rel.freq <- abs.freq/r
  ln.freq <- log(rel.freq)
  ln.freq[ln.freq == -Inf] <- 0
  
  #entropy
  entropy_matrix <- rel.freq * ln.freq
  entropy<- - colSums(entropy_matrix)
  norm_entropy <- entropy/log(num_cat)
  norm_entropy <- ifelse(is.nan(norm_entropy),0,norm_entropy)
  
  agreement <- vector(mode="numeric", length=s)
  ve <- matrix(data=0,nrow=r,ncol=r)
  row.names(ve) <- rnames
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        if (data[i,k] == data[j,k]) {
          agreement[k] <- norm_entropy[k]
        }
        else {
          agreement[k] <- 0
        }
      }
      if (i == j) {
        ve[i,j] <- 0
      }
      else {
        ve[i,j] <- 1-1/s*(sum(agreement))
        ve[j,i] <- ve[i,j]
      }
    }
  }
  return(as.dist(ve))
}