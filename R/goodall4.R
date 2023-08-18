#' Goodall 4 (G4) Measure
#' 
#' @description The function calculates a dissimilarity matrix based on the G4 similarity measure.
#'  
#' @param data A data.frame or a matrix with cases in rows and variables in columns.
#' 
#' @param var.weights A numeric vector setting weights to the used variables. One can choose the real numbers from zero to one.
#' 
#' @return The function returns an object of the class "dist".
#' \cr
#' 
#' @details The Goodall 4 similarity measure was presented in (Boriah et al., 2008). It is a simple modification of the original Goodall measure (Goodall, 1966).                                  
#' It assigns higher weights to the frequent categories matches.
#' 
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'  \cr
#'  \cr
#' Goodall V.D. (1966). A new similarity index based on probability. Biometrics, 22(4), p. 882.
#'
#' @seealso
#' \code{\link[nomclust]{anderberg}},
#' \code{\link[nomclust]{burnaby}},
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{gambaryan}},
#' \code{\link[nomclust]{goodall1}},
#' \code{\link[nomclust]{goodall2}},
#' \code{\link[nomclust]{goodall3}},
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
#' prox.goodall4 <- goodall4(data20)
#'
#' # dissimilarity matrix calculation with variable weights
#' weights.goodall4 <- goodall4(data20, var.weights = c(0.7, 1, 0.9, 0.5, 0))
#'
#' @export

goodall4 <- function(data, var.weights = NULL) {
  
  # dealing with the missing data
  if (sum(is.na(data)) > 0) {
    stop("The dissimilarity matrix CANNOT be calculated if the 'data' argument contains NA values.")
  }
  
  rnames <- row.names(data)
  
  # recoding everything to factors and then to numeric values
  indx <- sapply(data, is.factor)
  data[!indx] <- sapply(data[!indx], function(x) as.factor(x))
  data <- as.data.frame(unclass(data))
  data <- sapply(data, function(x) as.numeric(x))
  data <- as.data.frame(data)
  
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
  
  prox_matrix <- SIMILARITY(data, measure = "goodall4", freq.table, wt = var.weights)
  
  row.names(prox_matrix) <- rnames
  
  return(as.dist(prox_matrix))
}