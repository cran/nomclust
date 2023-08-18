#' Variable Entropy (VE) Measure
#' 
#' @description The function calculates a dissimilarity matrix based on the VE similarity measure.
#'  
#' @param data A data.frame or a matrix with cases in rows and variables in columns.
#' 
#' @param var.weights A numeric vector setting weights to the used variables. One can choose the real numbers from zero to one.
#' 
#' @return The function returns an object of the class "dist".
#' \cr
#' 
#' @details The Variable Entropy similarity measure was introduced in (Sulc and Rezankova, 2019). It treats
#' the similarity between two categories based on the within-cluster variability expressed by the normalized entropy.
#' The measure assigns higher weights to rare categories.
#' 
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'  \cr
#'  \cr
#' Sulc Z. and Rezankova H. (2019). Comparison of Similarity Measures for Categorical Data in Hierarchical Clustering. Journal of Classification. 2019, 35(1), p. 58-72. DOI: 10.1007/s00357-019-09317-5.
#'  
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
#' \code{\link[nomclust]{sm}},
#' \code{\link[nomclust]{smirnov}},
#' \code{\link[nomclust]{vm}}.
#'
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' # sample data
#' data(data20)
#' 
#' # dissimilarity matrix calculation
#' prox.ve <- ve(data20)
#' 
#' # dissimilarity matrix calculation with variable weighting
#' prox.ve.2 <- ve(data20, var.weights = c(1, 0.8, 0.6, 0.4, 0.2))
#' 
#' # dissimilarity matrix calculation with variable weights
#' weights.ve <- ve(data20, var.weights = c(0.7, 1, 0.9, 0.5, 0))
#' 
#' @export 

ve <- function(data, var.weights = NULL) {
  
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
  } else if (!(range(var.weights)[1] >= 0 & range(var.weights)[2] <= 1)) {
    stop("The weight vector should contain values in a range from zero to one.")
  }
  
  
  freq.table <- freq.abs(data)
  
  prox_matrix <- SIMILARITY(data, measure = "ve", freq.table, wt = var.weights)
  
  row.names(prox_matrix) <- rnames
  
  return(as.dist(prox_matrix))
}