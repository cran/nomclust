#' Cluster Quality Evaluation of Nominal Data Hierarchical Clustering
#' 
#' @description The function calculates a set of evaluation criteria if the original dataset and the cluster membership variables are provided. 
#'  The function calculates up to eight evaluation criteria described in (Sulc et al., 2018) and provides the optimal number of clusters based on these criteria. 
#'  It is primarily focused on evaluating hierarchical clustering results obtained by similarity measures different from those that occur in the nomclust package. 
#'  Thus, it can serve for the comparison of various similarity measures for categorical data.
#' 
#' @param data A data.frame or a matrix with cases in rows and variables in colums.
#' 
#' @param clusters A data.frame or a list of cluster memberships obtained based on the dataset defined in the parameter \code{data} in the form of a sequence from the two-cluster solution to the maximal-cluster solution.
#' 
#' @param diss An optional parameter. A matrix or a dist object containing dissimilarities calculated based on the dataset defined in the parameter \code{data}.
#' 
#' @return The function returns a list with three components.
#' \cr
#' \cr
#' The \code{eval} component contains up to eight evaluation criteria as vectors in a list. Namely, Within-cluster mutability coefficient (WCM), Within-cluster entropy coefficient (WCE),
#' Pseudo F Indices based on the mutability (PSFM) and the entropy (PSFE), Bayesian (BIC), and Akaike (AIC) information criteria for categorical data, the BK index, and, if the \code{diss.matrix} argument is present, the silhouette index (SI).
#' \cr
#' \cr
#' The \code{opt} component is present in the output together with the \code{eval} component. It displays the optimal number of clusters for the evaluation criteria from the \code{eval} component, except for WCM and WCE, where the optimal number of clusters is based on the elbow method.
#' \cr
#' \cr
#' The \code{call} component contains the function call.
#' 
#'@references
#' Sulc Z., Cibulkova J., Prochazka J., Rezankova H. (2018). Internal Evaluation Criteria for Categorical Data in Hierarchical Clustering: Optimal Number of Clusters Determination, Metodoloski Zveski, 15(2), p. 1-20.
#'
#' @seealso
#' \code{\link[nomclust]{nomclust}}, \code{\link[nomclust]{nomprox}}, \code{\link[nomclust]{eval.plot}}.
#' 
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' # sample data
#' data(data20)
#' 
#' # creating an object with results of hierarchical clustering
#' hca.object <- nomclust(data20, measure = "iof", method = "average", clu.high = 7)
#' 
#' # the cluster memberships
#' data20.clu <- hca.object$mem
#' 
#' # obtaining evaluation criteria for the provided dataset and cluster memberships
#' data20.eval <- evalclust(data20, clusters = data20.clu)
#' 
#' # visualization of the evaluation criteria
#' eval.plot(data20.eval)
#' 
#' # silhouette index can be calculated if the dissimilarity matrix is provided
#' data20.eval <- evalclust(data20, clusters = data20.clu, diss = hca.object$prox)
#' eval.plot(data20.eval, criteria = "SI")
#' 
#' @export 

evalclust <- function (data, clusters, diss = NULL) {
  
  clu_low = 2
  
  # change a possible list input to data.frame
  if (is.data.frame(clusters) == FALSE) {
    clusters <- as.data.frame(clusters)
  }
  
  # check the lenght of data
  if (nrow(data) != nrow(clusters)) {
    stop("The dataset and the cluster membership variables are of different lengths.")
  }
  
  if (ncol(clusters) < 2) {
    stop("The clusters argument must contain clusters for at least two cluster membership variables.")
  }
  
  # dealing with the missing data
  if (sum(is.na(data)) > 0) {
    stop("The evaluation CANNOT be run if the 'data' argument contains NA values.")
  }
  
  
  #if matrix, coerce to data.frame
  if(is.matrix(data) == 1) {
    data <- as.data.frame(data)
  }
  
  
  #number of variables of dataset
  num_var <- ncol(data)
  
  #max number of categories
  num_cat <- sapply(data, function(x) length(unique(x)))
  max_num_cat <- max(num_cat)
  
  
  clu.high <- ncol(clusters) + 1
  
  num_clu <- sapply(clusters, function(x) length(unique(x)))
  
  # to check the clusters object
  if (min(num_clu != 2)) {
    stop("The minimal number of clusters must be set to two.")
  }
  
  if (sum((num_clu - seq(clu_low, clu.high))^2) != 0) {
    clusters <- clusters[ ,order(num_clu)]
    num_clu <- sapply(clusters, function(x) length(unique(x)))
    if (sum((num_clu - seq(clu_low, clu.high))^2) != 0) {
      stop("In the 'clusters' parameter, a sequence of cluster memberships is not provided.")
    }
  }
  
  
    #creation of set of 3D matrices
    M <- list()
    for (i in clu_low:clu.high) {
      A <- list()
      A1 <- list()
      MMM <- array(0,dim=c(max_num_cat,i,num_var))
      M1 <- array(0,dim=c(max_num_cat,1,num_var))
      
      for (j in 1:num_var) {
        A[[j]] <- table(data[, j], clusters[,i - clu_low + 1])
        A1[[j]] <- rowSums(A[[j]])
      }
      
      for (j in 1:num_var) {
        MMM[1:nrow(A[[j]]), 1:ncol(A[[j]]), j] <- A[[j]]
        M1[1:nrow(A[[j]]),,j] <- A1[[j]]
      }
      M[[i-clu_low+2]] <- MMM
    }
    
    # to include one-cluster solution into the matrix
    M[[1]] <- M1
    
    #evaluation results
    results <- EVAL(M, clusters, diss)
    results1 <- results[[1]]
    results2 <- results[[2]]
    
    call <- match.call()
    
  list <- list(eval = results1, opt = results2, call = call)
  attr(list,"class")="nomclust"
  
  return(list)
}
