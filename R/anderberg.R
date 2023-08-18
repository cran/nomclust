#' Anderberg (AN) Measure
#' 
#' @description The function calculates a dissimilarity matrix based on the AN similarity measure.
#' 
#'  
#' @param data A data.frame or a matrix with cases in rows and variables in columns.
#' 
#' 
#' @return The function returns an object of the class "dist".
#' \cr
#' 
#' @details The Anderberg similarity measure was presented in (Anderberg, 1973).
#' The measure assigns higher weights to infrequent matches and mismatches. 
#' It takes on values from zero to one. The minimum similarity is attained when there are no matches and vice versa, see (Borian et al., 2008).
#' 
#' @references
#' Andergerg M.R. (1973). Cluster analysis for applications. Academic Press, New York.
#'  \cr
#'  \cr
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'
#' @seealso
#' \code{\link[nomclust]{burnaby}},
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{gambaryan}},
#' \code{\link[nomclust]{goodall1}},
#' \code{\link[nomclust]{goodall2}},
#' \code{\link[nomclust]{goodall3}},
#' \code{\link[nomclust]{goodall4}},
#' \code{\link[nomclust]{iof}},
#' \code{\link[nomclust]{lin}},
#' \code{\link[nomclust]{lin1}},
#' \code{\link[nomclust]{of}},
#' \code{\link[nomclust]{sm}},
#' \code{\link[nomclust]{smirnov}},
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
#' prox.anderberg <- anderberg(data20)
#' 
#' @export 

anderberg <- function(data) {
  
  # dealing with the missing data
  if (sum(is.na(data)) > 0) {
    stop("The dissimilarity matrix CANNOT be calculated if the 'data' argument contains NA values.")
  }
  
  rnames <- row.names(data)
  
  # recoding everything to factors and then to numeric values
  indx <- sapply(data, is.factor)
  data[!indx] <- lapply(data[!indx], function(x) as.factor(x))
  data <- as.data.frame(sapply(data, function(x) as.numeric(x)))
  
  freq.table <- freq.abs(data)
  var.weights <- rep(1, ncol(data))
  
  prox_matrix <- SIMILARITY(data, measure = "anderberg", freq.table, wt = var.weights)
  
  row.names(prox_matrix) <- rnames
  
  return(as.dist(prox_matrix))
}