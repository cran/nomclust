#' Morlini and Zani's (MZ) Measure 
#' 
#' @description A function for calculation of a proximity (dissimilarity) matrix based on the MZ similarity measure.
#' \cr
#' 
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#' 
#' @return The function returns a dissimilarity matrix of the size \code{n x n}, where \code{n} is the number of objects in the original dataset in the argument \code{data}.
#' \cr
#'
#' @details The MZ measure was originally introduced by Morlini and Zani (2012) under the name S2. The S2 measure was proposed. It is based on a binary-transformed dataset, so the \bold{morlini} function must first create dummy-coded variables.
#' The measure uses relative frequencies of categories of binary-coded variables, and it assigns higher weights to infrequent categories.
#'
#' @references
#' Morlini I., Zani S. (2012). A new class of weighted similarity indices using polytomous variables. Journal of Classification, 29(2), p. 199-226.
#' 
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
#' prox.morlini <- morlini(data20)
#' 
#' @export 




morlini <- function(data) {
  
  s <- ncol(data)
  #num_cat <- sapply(data, function(x) length(unique(x)))
  
  rnames <- row.names(data)
  
  # recoding everything to factors and then to numeric values
  indx <- sapply(data, is.factor)
  data[!indx] <- sapply(data[!indx], function(x) as.factor(x))
  data <- as.data.frame(unclass(data))
  data <- sapply(data, function(x) as.numeric(x))
  data <- as.data.frame(data)
  
  #with dummies
  #data_dummy <- dummy.data.frame(data, dummy.classes ="ALL",omit.constants = F)
  
  # dummy transformation of a data.frame
  data_dummy <- data.frame(row.names = row.names(data))
  num_cat <- sapply(data, function(x) length(unique(x)))
  for (i in 1:length(num_cat)) {
    variable_set <- data.frame()
    for (j in 1:num_cat[i]) {
      for (k in (1:nrow(data))) {
        variable_set[k,j] <- ifelse(data[k,i] == j, 1, 0)
      }
    }
    data_dummy <- cbind(data_dummy, variable_set)
  }
  
  n <- nrow(data_dummy)
  hs <- ncol(data_dummy)
  
  nsv <- sapply(data_dummy, sum)
  fsv2 <- log(1/(nsv/n)^2)
  
  E <- matrix(data=0,nrow=n,ncol=n)
  agreement <- vector(mode="numeric", length=hs)
  
  #computation of Eij
  for (i in 1:(n-1)) {
    for (j in (1+i):n) {
      for (k in 1:hs) {
        if (data_dummy[i,k] == 1 & data_dummy[j,k] == 1) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- 0
        }
      }
      E[i,j] <- fsv2 %*% agreement
      E[j,i] <- E[i,j]
    }
  }

#computation of Fij
  cum <- cumsum(num_cat)
  F <- matrix(data=0,nrow=n,ncol=n)
  row.names(F) <- rnames
  
  for (i in 1:(n-1)) {
    for (j in (1+i):n) {
      v <- 0
      agreement <- vector(mode="numeric", length=hs)
      for (k in 1:s) {
        for (t in (v+1):cum[k]) {
          if (data_dummy[i,t] == 0 & data_dummy[j,t] == 1) {
            agreement[(v+1):cum[k]] <- 1
          }
        }
        v <- cum[k]      
      }
      F[i,j] <- t(agreement) %*% fsv2
      F[j,i] <- F[i,j]
    }
  }
   morlini <- 1 - E/(E+F)
   morlini <- ifelse(is.nan(morlini),0,morlini)
  return(morlini)
}
  


  
  