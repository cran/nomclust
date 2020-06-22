#' Hierarchical Cluster Analysis for Nominal Data Based on a Proximity Matrix
#' 
#' @description The \code{nomprox()} function performs hierarchical cluster analysis in situations when the proximity (dissimilarity) matrix was calculated externally. For instance, in a different R package, in an own-created function, or in other software.
#' It offers three linkage methods that can be used for categorical data. The obtained clusters can be evaluated by seven evaluation indices, see (Sulc et al., 2018).
#' 
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#' 
#' @param diss A proximity matrix or a dist object calculated from the dataset defined in a parameter \code{data}.
#' 
#' @param clu.high A \emph{numeric} value expressing the maximal number of cluster for which the cluster memberships variables are produced.
#' 
#' @param eval A \emph{logical} operator; if TRUE, evaluation of clustering results is performed.
#' 
#' @param method A \emph{character} string defining the clustering method. The following methods can be used: \code{"average"}, \code{"complete"}, \code{"single"}.
#' 
#' @return The function returns a list with up to three components:
#' \cr
#' \cr
#' The \code{mem} component contains cluster membership partitions for the selected numbers of clusters in the form of a \emph{list}.
#' \cr
#' \cr
#' The \code{eval} component contains seven evaluation criteria in as vectors in a \emph{list}. Namely, Within-cluster mutability coefficient (WCM), Within-cluster entropy coefficient (WCE),
#' Pseudo F Indices based on the mutability (PSFM) and the entropy (PSFE), Bayessian (BIC) and Akaike (AIC) information criteria for categorical data and the BK index.
#' To see them all in once, the form of a \emph{data.frame} is more appropriate.
#' \cr
#' \cr
#' The \code{opt} component is present in the output together with the \code{eval} component. It displays the optimal number of clusters for the evaluation criteria from the \code{eval} component, except for WCM and WCE, where the optimal number of clusters is based on the elbow method.
#' 
#' @references
#' Sulc Z., Cibulkova J., Prochazka J., Rezankova H. (2018). Internal Evaluation Criteria for Categorical Data in Hierarchical Clustering: Optimal Number of Clusters Determination, Metodoloski Zveski, 15(2), p. 1-20.
#' 
#' @seealso
#' \code{\link[nomclust]{nomclust}}, \code{\link[nomclust]{evalclust}}, \code{\link[nomclust]{eval.plot}}.
#' 
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' # sample data
#' data(data20)
#' 
#' # computation of a dissimilarity matrix using the iof similarity measure
#' diss.matrix <- iof(data20)
#' 
#' # creating an object with results of hierarchical clustering 
#' hca.object <- nomprox(diss = diss.matrix, data = data20, method = "complete",
#'  clu.high = 5, eval = TRUE)
#' 
#' 
#' @export


nomprox <- function (diss, data = NULL, method = "average", clu.high = 6, eval = TRUE) {
  
  clu.low = 2
  
  # elimination of other clustering methods
  if (method %in% c("single", "average", "complete") == FALSE) {
    stop("An invalid clustering method was chosen.")
  }
  
  if (is.null(data) == TRUE & eval == TRUE) {
    eval <- FALSE
    warning("The argument 'eval' was set to FALSE since the 'data' argument needed for evaluation criteria calculation was not provided.")
  }
  
  #number of clusters cannot exceed the parameter clu.high
  if (nrow(diss)<clu.high) {
    stop("The argument 'clu.high' cannot exceed the number of clustered objects.")
  }
  
  # transforms the dist object into a matrix
  if (is(diss, "dist") == TRUE) {
    diss <- as.matrix(diss)
  }
  
  # is an argument a square proximity matrix
  if ((nrow(diss) != ncol(diss)) == TRUE) {
        stop("The argument 'diss' is not a square proximity matrix.")
  }
  
  # dealing with the missing data
  if (sum(is.na(diss)) > 0) {
    stop("The dissimilarity matrix contains NA values. It is probably damaged.")
  }
  
  
  if (eval == 1) {
    
    # check if the data dimensions correspond to the proximity matrx dimensions
    if (nrow(data) != nrow(diss)) {
      stop("The used dataset and the dissimilarity matrix are of different sizes.")
    }
    
    # dealing with the missing data
    if (sum(is.na(data)) > 0) {
      stop("The cluster analysis CANNOT be run if the 'data' argument contains NA values.")
    }
    
    
    # taking row.names from data
    rnames <- row.names(data)
    
    #if matrix, coerce to data.frame
    if(is.matrix(data) == 1) {
      data <- as.data.frame(data)
    }
    
    # recoding everything to factors and then to numeric values
    indx <- sapply(data, is.factor)
    data[!indx] <- sapply(data[!indx], function(x) as.factor(x))
    data <- as.data.frame(unclass(data))
    data <- sapply(data, function(x) as.numeric(x))
    data <- as.data.frame(data)
    
    #number of variables of dataset
    num_var <- ncol(data)
    
    #max number of categories
    num_cat <- sapply(data, function(x) length(unique(x)))
    max_num_cat <- max(num_cat)
    
    # frequencies of categories in all variables
    freq.table <- freq.abs(data)
    
    # adding row.names to proximity matrix
    row.names(diss) <- rnames
  }
  
  #hierarchical cluster analysis, where "prox" is a proximity matrix
  hca <- agnes(diss, diss = TRUE, method = method)
  #hca2 <- agnes(aaa$prox, diss = TRUE, method = method)
  
  #cluster membership
  data_clu <- vector(mode="numeric", length = nrow(diss))
  clu_name <- vector(mode="character", length = length(clu.low:clu.high)+1)
  for (i in clu.low:clu.high) {
    clusters <- cutree(hca, i)
    data_clu <-data.frame(data_clu, clusters)
    clu_name[i] <- paste("clu_", i, sep = "" )
  }
  names(data_clu) <- clu_name
  clusters <- data_clu[,-1]
  
  if (eval == 1) {
    #creation of set of 3D matrices
    M <- list()
    for (i in clu.low:clu.high) {
      A <- list()
      A1 <- list()
      MMM <- array(0,dim=c(max_num_cat,i,num_var))
      M1 <- array(0,dim=c(max_num_cat,1,num_var))
      
      for (j in 1:num_var) {
        A[[j]] <- table(data[, j], clusters[,i - clu.low + 1])
        A1[[j]] <- rowSums(A[[j]])
      }
      
      for (j in 1:num_var) {
        MMM[1:nrow(A[[j]]), 1:ncol(A[[j]]), j] <- A[[j]]
        M1[1:nrow(A[[j]]),,j] <- A1[[j]]
      }
      M[[i-clu.low+2]] <- MMM
    }
    
    # to include one-cluster solution into the matrix
    M[[1]] <- M1
    
    #evaluation results
    results <- EVAL(M)
    results1 <- results[[1]]
    results2 <- results[[2]]
    
  }
  
  clu_results <-  as.list(clusters)
  dend <- hca[-c(5,6)]
  
  if (eval == 1) {
    list <- list(mem = clu_results, eval = results1, opt = results2, dend = dend)
  }
  if (eval == 0) {
    list <- list(mem = clu_results, dend = dend)
  }
  
  
  return(list)
}
