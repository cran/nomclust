#' Goodall 3 (G3) Measure
#' 
#' @description The function calculates a dissimilarity matrix based on the G3 similarity measure.
#'  
#' @param data A data.frame or a matrix with cases in rows and variables in colums.
#' 
#' @return The function returns an object of the class "dist".
#' \cr
#' 
#' @details The Goodall 3 similarity measure was presented in (Boriah et al., 2008). It is a simple modification of the original Goodall measure (Goodall, 1966).           
#' The measure assigns higher weight if the infrequent categories match
#' regardless on frequencies of other categories.
#' 
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#' \cr
#' \cr
#' Goodall V.D. (1966). A new similarity index based on probability. Biometrics, 22(4), p. 882.
#'
#' @seealso
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{good1}},
#' \code{\link[nomclust]{good2}},
#' \code{\link[nomclust]{good4}},
#' \code{\link[nomclust]{iof}},
#' \code{\link[nomclust]{lin}},
#' \code{\link[nomclust]{lin1}},
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
#' prox.good3 <- good3(data20)
#'
#' @export 

good3 <- function(data) {
  
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
  # } else if (is.numeric(var.weights) == TRUE) {
  #    if(is.na(sum(var.weights >= 0)) | sum(var.weights >= 0)!=ncol(data)) {
  #     stop("The vector of weights contains negative or missing values.")
  #  }
  #    var.wgt <- var.weights
  
  
  # } else {
  #   stop("Invalid weighting scheme.")
  # }
  freq.table <- freq.abs(data)

  prox_matrix <- SIMILARITY(data, measure = "good3", freq.table)
  
  row.names(prox_matrix) <- rnames
  
  return(as.dist(prox_matrix))
}