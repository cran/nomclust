#' Inverse Occurence Frequency (IOF) Measure
#' 
#' @description A function for calculation of a proximity (dissimilarity) matrix based on the IOF similarity measure.
#' \cr               
#'  
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#'
#' @return The function returns a dissimilarity matrix of the size \code{n x n}, where \code{n} is the number of objects in the original dataset in the argument \code{data}.
#' \cr
#'
#' @details The IOF (Inverse Occurrence Frequency) measure was originally constructed for the text mining tasks,
#' see (Sparck-Jones, 1972), later, it was adjusted for categorical variables, see (Boriah et al., 2008).
#' The measure assigns higher weight to mismatches on less frequent values and vice versa. 
#'
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'  \cr
#'  \cr
#' Spark-Jones K. (1972). A statistical interpretation of term specificity and its application in retrieval.
#' In Journal of Documentation, 28(1), 11-21. Later: Journal of Documentation, 60(5) (2002), 493-502.
#' 
#' @seealso
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{good1}},
#' \code{\link[nomclust]{good2}},
#' \code{\link[nomclust]{good3}},
#' \code{\link[nomclust]{good4}},
#' \code{\link[nomclust]{lin}},
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
#' prox.iof <- iof(data20)
#' 
#' @export 



iof <- function(data) {
  
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
  
  agreement <- vector(mode="numeric", length=s)
  iof <- matrix(data=0,nrow=r,ncol=r)
  row.names(iof) <- rnames
  
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        d <- data[j,k]
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- 1/(1+log(freq.abs[c,k])*log(freq.abs[d,k]))
        }
      }
      iof[i,j] <- 1/(1/s*(sum(agreement)))-1
      iof[j,i] <- iof[i,j]
    }
  }
  return(iof)
}