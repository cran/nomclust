#' Occurence Frequency (OF) Measure
#' 
#' @description A function for calculation of a proximity (dissimilarity) matrix based on the OF similarity measure.
#'                                        
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#' 
#' @return The function returns an object of class "dist".
#' \cr
#'
#' @details The OF (Occurrence Frequency) measure was originally constructed for the text mining tasks,
#' see (Sparck-Jones, 1972), later, it was adjusted for categorical variables, see (Boriah et al., 2008)
#' It assigns higher weight to mismatches on less frequent values and otherwise.
#'
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'  \cr
#'  \cr
#' Spark-Jones K. (1972). A statistical interpretation of term specificity and its application in retrieval.
#' In Journal of Documentation, 28(1), p. 11-21. Later: Journal of Documentation, 60(5) (2002), p. 493-502.
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
#' prox.of <- of(data20)
#' 
#' @export 




of <- function(data) {
  
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
  

  freq.abs <- freq.abs(data)
  
  agreement <- vector(mode="numeric", length=s)
  of <- matrix(data=0,nrow=r,ncol=r)
  row.names(of) <- rnames
  
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        d <- data[j,k]
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- 1/(1+log(r/freq.abs[c,k])*log(r/freq.abs[d,k]))
        }
      }
      of[i,j] <- 1/(1/s*(sum(agreement)))-1
      of[j,i] <- of[i,j]
    }
  }
  return(as.dist(of))
}

