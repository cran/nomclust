#' Simple Matching Coefficient (SM)
#' 
#' @description A function for calculation of a proximity (dissimilarity) matrix based on the SM similarity measure.
#'  
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#' 
#' @return The function returns a dissimilarity matrix of the size \code{n x n}, where \code{n} is the number of objects in the original dataset in the argument \code{data}.
#' \cr
#' 
#' @details The simple matching coefficient (Sokal, 1958) represents the simplest way of measuring similarity. It does not impose any weights.
#' By a given variable, it assigns the value 1 in case of match and value 0 otherwise.
#' 
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'  \cr
#'  \cr
#' Sokal R., Michener C. (1958). A statistical method for evaluating systematic relationships. In: Science bulletin, 38(22),
#' The University of Kansas.  
#'  \cr
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
#' prox.sm <- sm(data20)
#'
#' @export 



sm <- function(data) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  rnames <- row.names(data)
  
  # recoding everything to factors and then to numeric values
  indx <- sapply(data, is.factor)
  data[!indx] <- sapply(data[!indx], function(x) as.factor(x))
  data <- as.data.frame(unclass(data))
  data <- sapply(data, function(x) as.numeric(x))
  data <- as.data.frame(data)
  
  agreement <- vector(mode="numeric", length=s)
  sm <- matrix(data=0,nrow=r,ncol=r)
  row.names(sm) <- rnames
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- 0
        }
      }
      sm[i,j] <- 1-1/s*(sum(agreement))
      sm[j,i] <- sm[i,j]
    }
  }
  return(sm)
}
