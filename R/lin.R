#' Lin (LIN) Measure
#' 
#' @description A function for calculation of a proximity (dissimilarity) matrix based on the LIN similarity measure.
#' 
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#' 
#' @return The function returns a dissimilarity matrix of the size \code{n x n}, where \code{n} is the number of objects in the original dataset in the argument \code{data}.
#' \cr
#'
#' @details The Lin measure was introduced by Lin (1998) and presented in (Boriah et al., 2008).
#' The measure assigns higher weights to more frequent categories in case of matches
#' and lower weights to less frequent categories in case of mismatches.
#'
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'  \cr
#'  \cr
#' Lin D. (1998). An information-theoretic definition of similarity.
#' In: ICML '98: Proceedings of the 15th International Conference on Machine Learning. San Francisco, p. 296-304.
#' 
#' @seealso
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{good1}},
#' \code{\link[nomclust]{good2}},
#' \code{\link[nomclust]{good3}},
#' \code{\link[nomclust]{good4}},
#' \code{\link[nomclust]{iof}},
#' \code{\link[nomclust]{lin1}},
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
#' prox.lin <- lin(data20)
#' 
#' @export 




lin <- function(data) {
  
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

  agreement <- vector(mode="numeric", length=s)
  lin <- matrix(data=0,nrow=r,ncol=r)
  row.names(lin) <- rnames
  weights <- vector(mode="numeric", length=s)
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        d <- data[j,k]
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 2*log(freq.rel[c,k])
        }
        else {
          agreement[k] <- 2*log(freq.rel[c,k] + freq.rel[d,k])
        }
        weights[k] <- log(freq.rel[c,k]) + log(freq.rel[d,k])
      }
      if (i == j) {
        lin[i,j] <- 0
      }
      else {
        lin[i,j] <- 1/(1/sum(weights)*(sum(agreement))) - 1
        lin[j,i] <- lin[i,j]
      }
    }
  }
  lin[lin == -Inf] <- max(lin) + 1
  return(lin)
}