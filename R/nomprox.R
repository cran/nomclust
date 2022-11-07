#' Hierarchical Clustering of Nominal Data Based on a Proximity Matrix
#' 
#' @description The function performs hierarchical cluster analysis based on a dissimilarity matrix.
#' 
#' @param data A data.frame or a matrix with cases in rows and variables in columns.
#' 
#' @param diss A proximity matrix or a dist object calculated based on the dataset defined in a parameter \code{data}.
#' 
#' @param clu.high A numeric value that expresses the maximal number of clusters for which the cluster membership variables are produced.
#' 
#' @param eval A logical operator; if TRUE, evaluation of clustering results is performed.
#' 
#' @param method A character string defining the clustering method. The following methods can be used: \code{"average"}, \code{"complete"}, \code{"single"}.
#' 
#' @param prox A logical operator or a numeric value. If a logical value TRUE indicates that the proximity matrix is a part of the output. A numeric value (integer) of this argument indicates the maximal number of cases in a dataset for which a proximity matrix will occur in the output.
#' 
#' @return The function returns a list with up to six components:
#' \cr
#' \cr
#' The \code{mem} component contains cluster membership partitions for the selected numbers of clusters in the form of a list.
#' \cr
#' \cr
#' The \code{eval} component contains up to eight evaluation criteria as vectors in a list. Namely, Within-cluster mutability coefficient (WCM), Within-cluster entropy coefficient (WCE),
#' Pseudo F Indices based on the mutability (PSFM) and the entropy (PSFE), Bayesian (BIC), and Akaike (AIC) information criteria for categorical data, the BK index, and, if the prox component is present, the silhouette index (SI).
#' \cr
#' \cr
#' The \code{opt} component is present in the output together with the \code{eval} component. It displays the optimal number of clusters for the evaluation criteria from the \code{eval} component, except for WCM and WCE, where the optimal number of clusters is based on the elbow method.
#' \cr
#' \cr
#' The \code{dend} component can be found in the output only together with the \code{prox} component. It contains all the necessary information for dendrogram creation.
#' \cr
#' \cr
#' The \code{prox} component contains the dissimilarity matrix in the form of the "dist" object.
#' \cr
#' \cr
#' The \code{call} component contains the function call.
#' 
#' 
#' @details The function performs hierarchical cluster analysis in situations when the proximity (dissimilarity) matrix was calculated externally. For instance, in a different R package, in an own-created function, or in other software.
#' It offers three linkage methods that can be used for categorical data. The obtained clusters can be evaluated by up to eight evaluation indices (Sulc et al., 2018).
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
#'  clu.high = 5, eval = TRUE, prox = FALSE)
#'  
#' # quick clustering summary
#' summary(hca.object)
#' 
#' # quick cluster quality evaluation
#' print(hca.object)
#' 
#' # visualization of the evaluation criteria
#' eval.plot(hca.object)
#' 
#' # a dendrogram can be displayed if the object contains the prox component
#' hca.object <- nomprox(diss = diss.matrix, data = data20, method = "complete",
#'  clu.high = 5, eval = TRUE, prox = TRUE)
#' 
#' # a quick dendrogram
#' plot(hca.object)
#' 
#' # a dendrogram with three designated clusters
#' dend.plot(hca.object, clusters = 3)
#' 
#' 
#' @export


nomprox <- function (diss, data = NULL, method = "average", clu.high = 6, eval = TRUE, prox = 100) {
  
  clu.low = 2
  
  # elimination of other clustering methods
  if (method %in% c("single", "average", "complete") == FALSE) {
    stop("An invalid clustering method was chosen.")
  }
  
  if (is.null(data) == TRUE & eval == TRUE) {
    eval <- FALSE
    warning("The argument 'eval' was set to FALSE since the 'data' argument needed for evaluation criteria calculation was not provided.")
  }
  
  # transforms the dist object into a matrix
  if (is(diss, "dist") == TRUE) {
    diss <- as.matrix(diss)
  }
  
  #number of clusters cannot exceed the parameter clu.high
  if (nrow(diss)<clu.high) {
    stop("The argument 'clu.high' cannot exceed the number of clustered objects.")
  }
  
  if (clu.high < 3) {
    stop("The 'clu.high' argument cannot be lower than 3.")
  }
  
  # calculate proximity matrix for up to maximal size of a dataset
  if (prox != FALSE & prox != TRUE) {
    if (is.numeric(prox) == T) {
      if (nrow(diss) <= abs(as.integer(prox))) {
        prox <- T
      } else
        prox <- F
    } else
      stop("The 'prox' argument should be of the 'numeric' type.")
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
    
    # check if the data dimensions correspond to the proximity matrix dimensions
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
    data[!indx] <- lapply(data[!indx], function(x) as.factor(x))
    data <- as.data.frame(sapply(data, function(x) as.numeric(x)))
    
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
    results <- EVAL(M, clusters, diss = diss)
    results1 <- results[[1]]
    results2 <- results[[2]]
    
  }
  
  clu_results <-  as.list(clusters)
  dend <- hca[-c(5,6)]
  
  call <- match.call()
  
  if (prox == 1) {
    if (eval == 1) {
      list <- list(mem = clu_results, eval = results1, opt = results2, dend = dend, prox = as.dist(diss), call = call)
    }
    if (eval == 0) {
      list <- list(mem = clu_results, dend = dend, prox = as.dist(diss), call = call)
    }
  }
  
  if (prox == 0) {
    if (eval == 1) {
      list <- list(mem = clu_results, eval = results1, opt = results2, dend = dend, call = call)
    }
    if (eval == 0) {
      list <- list(mem = clu_results, dend = dend, call = call)
    }
  }
  
  
  attr(list,"class")="nomclust"
  
  return(list)
}
