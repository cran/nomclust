#' Eskin (ES) Measure
#' 
#' @description The function calculates a dissimilarity matrix based on the ES similarity measure.
#' 
#' @param data A data.frame or a matrix with cases in rows and variables in colums.
#' 
#' @return The function returns an object of the class "dist".
#' \cr
#' 
#' @details The Eskin similarity measure was proposed by Eskin et al. (2002) and examined by Boriah et al., (2008). It is constructed to assign
#' higher weights to mismatches on variables with more categories.
#' 
#' @references
#' Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#'  \cr
#'  \cr
#' Eskin E., Arnold A., Prerau M., Portnoy L. and Stolfo S. (2002). A geometric framework for unsupervised anomaly detection.
#' In D. Barbara and S. Jajodia (Eds): Applications of Data Mining in Computer Security, p. 78-100. Norwell: Kluwer Academic Publishers.
#' 
#' @seealso
#' \code{\link[nomclust]{good1}},
#' \code{\link[nomclust]{good2}},
#' \code{\link[nomclust]{good3}},
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
#' prox.eskin <- eskin(data20)
#' @export 

eskin <- function(data) {
  
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
  
  prox_matrix <- SIMILARITY(data, measure = "eskin")
  
  row.names(prox_matrix) <- rnames
  
  return(as.dist(prox_matrix))
}
