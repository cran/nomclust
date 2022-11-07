#' Lin 1 (LIN1) Measure
#' 
#' @description The function calculates a dissimilarity matrix based on the LIN1 similarity measure.
#'                                        
#' @param data A data.frame or a matrix with cases in rows and variables in columns.
#' 
#' @param var.weights A numeric vector setting weights to the used variables. One can choose the real numbers from zero to one.
#' 
#' @return The function returns an object of the class "dist".
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
#' \code{\link[nomclust]{anderberg}},
#' \code{\link[nomclust]{burnaby}},
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{gambaryan}},
#' \code{\link[nomclust]{good1}},
#' \code{\link[nomclust]{goodall1}},
#' \code{\link[nomclust]{good2}},
#' \code{\link[nomclust]{goodall2}},
#' \code{\link[nomclust]{good3}},
#' \code{\link[nomclust]{goodall3}},
#' \code{\link[nomclust]{good4}},
#' \code{\link[nomclust]{goodall4}},
#' \code{\link[nomclust]{iof}},
#' \code{\link[nomclust]{lin}},
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
#' prox.lin1 <- lin1(data20)
#'
#' # dissimilarity matrix calculation with variable weights
#' weights.lin1 <- lin1(data20, var.weights = c(0.7, 1, 0.9, 0.5, 0))
#'
#' @export 


lin1 <- function(data, var.weights = NULL) {
  
  # dealing with the missing data
  if (sum(is.na(data)) > 0) {
    stop("The dissimilarity matrix CANNOT be calculated if the 'data' argument contains NA values.")
  }
  
  rnames <- row.names(data)
  
  # recoding everything to factors and then to numeric values
  indx <- sapply(data, is.factor)
  data[!indx] <- lapply(data[!indx], function(x) as.factor(x))
  data <- as.data.frame(sapply(data, function(x) as.numeric(x)))
  
  # variable weighting
  
  # if (var.weights %in% c("none", "MI", "nMI", "MU", "MA") == TRUE) {
  #   var.wgt <- WGT(data, var.weights, alpha)
  
  # OWN-DEFINED WEIGHTS
  if (is.null(var.weights) == TRUE) {
    var.weights <- rep(1, ncol(data))
  } else if (!(is.numeric(var.weights) & length(var.weights) == ncol(data))) {
    stop("The weight vector should be numeric with the length equal to the number of clustered variables.")
  } else if (!all(is.finite(var.weights))) {
    stop("The weight vector can contain only finite numbers in a range from zero to one.")
  } else if (!(range(var.weights)[1] >= 0 & range(var.weights)[2] <= 1)) {
    stop("The weight vector should contain values in a range from zero to one.")
  }
  
  
  freq.table <- freq.abs(data)
  
  prox_matrix <- SIMILARITY(data, measure = "lin1", freq.table, wt = var.weights)
  
  row.names(prox_matrix) <- rnames
  
  return(as.dist(prox_matrix))
}