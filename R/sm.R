#' Simple Matching Coefficient (SM)
#' 
#' @description The function calculates a dissimilarity matrix based on the SM similarity measure.
#'  
#' @param data A data.frame or a matrix with cases in rows and variables in columns.
#' 
#' @param var.weights A numeric vector setting weights to the used variables. One can choose the real numbers from zero to one.
#' 
#' @return The function returns an object of the class "dist".
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
#' \code{\link[nomclust]{anderberg}},
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
#' prox.sm <- sm(data20)
#'
#' # dissimilarity matrix calculation with variable weights
#' weights.sm <- sm(data20, var.weights = c(0.7, 1, 0.9, 0.5, 0))
#'
#' @export 



sm <- function(data, var.weights = NULL) {
  
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
  
  
  prox_matrix <- SIMILARITY(data, measure = "sm", freq.table = NULL, wt = var.weights)
  
  row.names(prox_matrix) <- rnames
  
  return(as.dist(prox_matrix))
}
