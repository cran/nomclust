#' Lin 1 (LIN1) Measure
#' 
#' @description A function for calculation of a proximity (dissimilarity) matrix based on the LIN1 similarity measure.
#'                                        
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#' 
#' @return The function returns a dissimilarity matrix of the size \code{n x n}, where \code{n} is the number of objects in the original dataset in the argument \code{data}.
#' \cr
#'
#' @details The Lin 1 similarity measure was introduced in (Boriah et al., 2008) as a modification of the original Lin measure (Lin, 1998). In has
#' a complex system of weights. In case of mismatch, lower similarity is assigned if either
#' the mismatching values are very frequent or their relative frequency is in between the relative
#' frequencies of mismatching values. Higher similarity is assigned if the mismatched categories
#' are infrequent and there are a few other infrequent categories. In case of match,
#' lower similarity is given for matches on frequent categories or matches on categories
#' that have many other values of the same frequency. Higher similarity is given to matches
#' on infrequent categories.
#' 
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'  \cr
#'  \cr
#' Lin D. (1998). An information-theoretic definition of similarity.
#' In: ICML '98: Proceedings of the 15th International Conference on Machine Learning. San Francisco, p. 296-304.
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
#' \code{\link[nomclust]{morlini}},
#' \code{\link[nomclust]{of}},
#' \code{\link[nomclust]{sm}},
#' \code{\link[nomclust]{ve}},
#' \code{\link[nomclust]{vm}}.
#'
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' # sample data
#' data(data20)
#' 
#' # dissimilarity matrix calculation
#' prox.lin1 <- lin1(data20)
#'
#' @export 


lin1 <- function(data) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  rnames <- row.names(data)
  
  # recoding everything to factors and then to numeric values
  indx <- sapply(data, is.factor)
  data[!indx] <- sapply(data[!indx], function(x) as.factor(x))
  data <- as.data.frame(unclass(data))
  data <- sapply(data, function(x) as.numeric(x))
  data <- as.data.frame(data)
  
  
  freq.abs <- freq.abs(data)
  freq.rel <- freq.abs/r
  freq.ln <- log(freq.rel)
  freq.ln[freq.ln == -Inf] <- 0
  
  
  agreement <- vector(mode="numeric", length=s)
  lin1 <- matrix(data=0,nrow=r,ncol=r)
  row.names(lin1) <- rnames
  weights <- vector(mode="numeric", length=s)

  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        d <- data[j,k]
        if (data[i,k] == data[j,k]) {
          logic <- freq.rel[,k] == freq.rel[c,k]
          agreement[k] <- sum(logic * freq.ln[,k])
          weights[k] <- sum(logic * freq.ln[,k])
        }
        else {
          if (freq.rel[c,k] >= freq.rel[d,k]) {
            logic <- freq.rel[,k] >= freq.rel[d,k] & freq.rel[,k] <= freq.rel[c,k]
            agreement[k] <- 2*log(sum(logic * freq.rel[,k]))
            weights[k] <- sum(logic * freq.ln[,k])
          }
          else {
            logic <- freq.rel[,k] >= freq.rel[c,k] & freq.rel[,k] <= freq.rel[d,k]
            agreement[k] <- 2*log(sum(logic * freq.rel[,k]))
            weights[k] <- sum(logic * freq.ln[,k])
          }
        }
      }
      lin1[i,j] <- 1/(1/sum(weights)*(sum(agreement))) - 1
      lin1[j,i] <- lin1[i,j]
    }
  }
  lin1[lin1 == -Inf] <- max(lin1) + 1
  return(lin1)
}