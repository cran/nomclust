#' Gambaryan (GA) Measure
#' 
#' @description The function calculates a dissimilarity matrix based on the GA similarity measure.
#' 
#'  
#' @param data A data.frame or a matrix with cases in rows and variables in columns.
#' 
#' 
#' @return The function returns an object of the class "dist".
#' \cr
#' 
#' @details The Gambaryan similarity measure was presented in (Gambaryan, 1964).
#' The measure assigns low weight to matches where the matching value occurs in about half the dataset, i.e., in between being frequent and rare, see (Borian et al., 2008).
#' 
#' @references
#' Gambaryan P. (1964). A mathematical model of taxonomy. 
#' SSR, 17(12), 47-53.
#'  \cr
#'  \cr
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'
#' @seealso
#' \code{\link[nomclust]{anderberg}},
#' \code{\link[nomclust]{burnaby}},
#' \code{\link[nomclust]{eskin}},
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
#' prox.gambaryan <- gambaryan(data20)
#' 
#' @export 

gambaryan <- function(data) {
  
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
  # if (is.null(var.weights) == TRUE) {
  #   var.weights <- rep(1, ncol(data))
  # } else if (!(is.numeric(var.weights) & length(var.weights) == ncol(data))) {
  #   stop("The weight vector should be numeric with the length equal to the number of clustered variables.")
  # } else if (!all(is.finite(var.weights))) {
  #   stop("The weight vector can contain only finite numbers in a range from zero to one.")
  # } else if (!(range(var.weights)[1] >= 0 & range(var.weights)[2] <= 1)) {
  #   stop("The weight vector should contain values in a range from zero to one.")
  # }
  
  
  freq.table <- freq.abs(data)
  var.weights <- rep(1, ncol(data))
  
  prox_matrix <- SIMILARITY(data, measure = "gambaryan", freq.table, wt = var.weights)
  
  row.names(prox_matrix) <- rnames
  
  return(as.dist(prox_matrix))
}