#' Convert Objects to Class agnes, twins
#' 
#' @description Converts objects of the class "nomclust" to the class "agnes, twins".
#' 
#' @param x The "nomclust" object containing components "dend" and "prox".
#' 
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return The function returns an object of class "agnes, twins".
#' \cr
#' 
#' @seealso
#' \code{\link[cluster]{agnes}}, \code{\link[stats]{as.hclust}} and \code{\link[stats]{hclust}}.
#' 
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' # sample data
#' data(data20)
#'
#' # creating an object with results of hierarchical clustering of 
#' hca.object <- nomclust(data20, measure = "lin", method = "average",
#'  clu.high = 5, prox = TRUE)
#' 
#' # nomclust plot
#' plot(hca.object)
#' 
#' # obtaining the agnes, twins object
#' hca.object.agnes <- as.agnes(hca.object)
#' 
#' # agnes plot
#' plot(hca.object.agnes)
#' 
#' # obtaining the hclust object
#' hca.object.hclust <- as.hclust(hca.object)
#' 
#' # hclust plot
#' plot(hca.object.hclust)
#' 
#' @export

as.agnes <- function(x, ...) {
  
  if (is(x, "nomclust") == TRUE) {
    
    if ((is.null(x$prox) == FALSE) & (is.null(x$dend) == FALSE)) { # nomclust object and nomprox with data
      
      output <- list(
        order = x$dend$order,
        height = x$dend$height,
        ac = x$dend$ac,
        merge = x$dend$merge,
        diss = x$prox,
        call = x$call,
        method = x$dend$method,
        order.lab = x$dend$order.lab
      )
      attr(output,"class") = c("agnes", "twins")
      return(output)
      
    } else {
      stop("Invalid object type. The function can be applied only on the outputs of nomclust() or nomprox() functions with 'dend' and 'prox' components.")
    }
  } else {
    stop("Invalid object type. The function can be applied only on the outputs of nomclust() or nomprox() functions with 'dend' and 'prox' components.")
  }

}
