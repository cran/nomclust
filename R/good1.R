#' Goodall 1 (G1) Measure
#' 
#' @description A function for calculation of a proximity (dissimilarity) matrix based on the G1 similarity measure.
#' 
#'  
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#' 
#' @return The function returns a dissimilarity matrix of the size \code{n x n}, where \code{n} is the number of objects in the original dataset in the argument \code{data}.
#' \cr
#' 
#' @details The Goodall 1 similarity measure was presented in (Boriah et al., 2008).
#' It is a simple modification of the original Goodall measure (Goodall, 1966). The measure assigns higher weights to infrequent matches.
#' 
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'  \cr
#'  \cr
#' Goodall V.D. (1966). A new similarity index based on probability. Biometrics, 22(4), p. 882.
#'
#' @seealso
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{good2}},
#' \code{\link[nomclust]{good3}},
#' \code{\link[nomclust]{good4}},
#' \code{\link[nomclust]{iof}},
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
#' prox.good1 <- good1(data20)
#' 
#' @export 

good1 <- function(data) {
  
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
  freq.rel2 <- freq.rel^2
  
  agreement <- vector(mode="numeric", length=s)
  good1 <- matrix(data=0,nrow=r,ncol=r)
  row.names(good1) <- rnames
  
  
  for (i in 1:(r-1)) {
    for (j in (1+i):r) {
      for (k in 1:s) {
        c <- data[i,k]
        if (data[i,k] == data[j,k]) {
          logic <- t(freq.rel[,k] <= freq.rel[c,k])
          agreement[k] <- 1 - sum(freq.rel2[,k] * logic)
        }
        else {
          agreement[k] <- 0
        }
      }
      if (i == j) {
        good1[i,j] <- 0
      }
      else {
        good1[i,j] <- 1-1/s*(sum(agreement))
        good1[j,i] <- good1[i,j]
      }
    }
  }
  return(good1)
}