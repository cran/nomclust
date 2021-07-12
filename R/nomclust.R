#' Hierarchical Clustering of Nominal Data
#' 
#' @description The function runs hierarchical cluster analysis (HCA) with objects characterized by nominal variables (without natural order of categories).
#'  It completely covers the clustering process, from the dissimilarity matrix calculation to the cluster quality evaluation. The function enables a user to choose from twelve similarity measures for nominal data summarized by (Boriah et al., 2008) and by (Sulc and Rezankova, 2019). 
#'  Next, it offers to choose from three linkage methods that can be used for categorical data. The obtained clusters can be evaluated by up to eight evaluation criteria (Sulc et al., 2018). The output of the nomclust() function may serve as an input for the visualization functions \emph{dend.plot} and \emph{eval.plot} in the nomclust package.
#' 
#' 
#' @param data A data.frame or a matrix with cases in rows and variables in columns.
#' 
#' @param measure A character string defining the similarity measure used for computation of proximity matrix in HCA:
#' \code{"eskin"}, \code{"good1"}, \code{"good2"}, \code{"good3"}, \code{"good4"}, \code{"iof"}, \code{"lin"}, \code{"lin1"}, \code{"of"}, \code{"sm"}, \code{"ve"}, \code{"vm"}.
#' 
#' @param method A character string defining the clustering method. The following methods can be used: \code{"average"}, \code{"complete"}, \code{"single"}.
#' 
#' @param clu.high  A numeric value expressing the maximal number of cluster for which the cluster memberships variables are produced.
#' 
#' @param eval A logical operator; if TRUE, evaluation of the clustering results is performed.
#' 
#' @param prox A logical operator or a numeric value. If a logical value TRUE indicates that the proximity matrix is a part of the output. A numeric value (integer) of this argument indicates the maximal number of cases in a dataset for which a proximity matrix will occur in the output.
#' 
#' 
#' @return The function returns a list with up to six components.
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
#' The \code{dend} component can be found in the output together with the \code{prox} component. It contains all the necessary information for dendrogram creation.
#' \cr
#' \cr
#' The \code{prox} component contains the dissimilarity matrix in the form of the "dist" object.
#' \cr
#' \cr
#' The \code{call} component contains the function call.
#' 
#'@references
#' Boriah S., Chandola V. and Kumar, V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#' \cr
#' \cr
#' Sulc Z., Cibulkova J., Prochazka J., Rezankova H. (2018). Internal Evaluation Criteria for Categorical Data in Hierarchical Clustering: Optimal Number of Clusters Determination, Metodoloski Zveski, 15(2), p. 1-20.
#' \cr
#' \cr
#' Sulc Z. and Rezankova H. (2019). Comparison of Similarity Measures for Categorical Data in Hierarchical Clustering. Journal of Classification. 2019, 35(1), p. 58-72. DOI: 10.1007/s00357-019-09317-5.
#' 
#' 
#' @seealso
#' \code{\link[nomclust]{evalclust}}, \code{\link[nomclust]{nomprox}}, \code{\link[nomclust]{eval.plot}}, \code{\link[nomclust]{dend.plot}}.
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
#' # quick clustering summary
#' summary(hca.object)
#' 
#' # quick cluster quality evaluation
#' print(hca.object)
#' 
#' # visualization of the evaluation criteria
#' eval.plot(hca.object)
#' 
#' # a quick dendrogram
#' plot(hca.object)
#' 
#' # a dendrogram with three designated clusters
#' dend.plot(hca.object, clusters = 3)
#' 
#' # obtaining values of evaluation indices as a data.frame
#' data20.eval <- as.data.frame(hca.object$eval)
#' 
#' # getting the optimal numbers of clusters as a data.frame
#' data20.opt <- as.data.frame(hca.object$opt)
#' 
#' # extracting cluster membership variables as a data.frame
#' data20.mem <- as.data.frame(hca.object$mem)
#' 
#' # obtaining a proximity matrix
#' data20.prox <- as.matrix(hca.object$prox)
#' 
#' # setting the maximal number of objects for which a proximity matrix is provided in the output to 30
#' hca.object <- nomclust(data20, measure = "iof", method = "complete",
#'  clu.high = 5, prox = 30)
#'  
#' # transforming the nomclust object to the class "hclust"
#' hca.object.hclust <- as.hclust(hca.object)
#' 
#' # transforming the nomclust object to the class "agnes, twins"
#' hca.object.agnes <- as.agnes(hca.object)
#' 
#' 
#' @export

nomclust <- function (data, measure = "lin", method = "average", clu.high = 6, eval = TRUE, 
                       prox = 100) {
  
  clu.low = 2
  
  # elimination of other clustering methods
  if (method %in% c("single", "average", "complete") == FALSE) {
    stop("An invalid clustering method was chosen.")
  }
  
  # legacy - morlini
  if (measure == "morlini") {
    stop("The 'morlini' similarity measure was removed from the package in version 2.2.1. 
       If you want to send a dissimilarity matrix calculation script using the 'morlini' measure, please contact the package's maintainer.")
  }
  
  # check of the used similarity measure
  if (measure %in% c("eskin", "good1", "good2", "good3", "good4", "iof", "of", "lin", "lin1", "sm", "ve", "vm") == FALSE) {
    stop("Invalid name of the similarity measure.")
  }
  
  # number of clusters cannot exceed the parameter clu.high
  if (nrow(data)<clu.high) {
    stop("The 'clu.high' argument cannot exceed the number of clustered objects.")
  }
  
  if (clu.high < 3) {
    stop("The 'clu.high' argument cannot be lower than 3.")
  }
  
  # calculate proximity matrix for up to maximal size of a dataset
  if (prox != FALSE & prox != TRUE) {
    if (is.numeric(prox) == T) {
      if (nrow(data) <= abs(as.integer(prox))) {
        prox <- T
      } else
        prox <- F
    } else
      stop("The 'prox' argument should be of the 'numeric' type.")
  }
  
  # dealing with the missing data
  if (sum(is.na(data)) > 0) {
    stop("The cluster analysis CANNOT be run if the 'data' argument contains NA values.")
  }
  
  rnames <- row.names(data)
  
  # if matrix, coerce to data.frame
  if(is.matrix(data) == 1) {
    data <- as.data.frame(data)
  }
  
  # recoding everything to factors and then to numeric values
  indx <- sapply(data, is.factor)
  data[!indx] <- lapply(data[!indx], function(x) as.factor(x))
  data <- as.data.frame(sapply(data, function(x) as.numeric(x)))
  
  
  # number of variables of dataset
  num_var <- ncol(data)
  
  # frequency distribution table
  freq.table <- freq.abs(data)
  
  # max number of categories
  num_cat <- sapply(data, function(x) length(unique(x)))
  max_num_cat <- max(num_cat)
  
  #computing the proximity matrices
  diss.matrix <- SIMILARITY(data, measure, freq.table) # casem pridat var.weight argument
  row.names(diss.matrix) <- rnames
  
  #hierarchical cluster analysis, where "prox" is a proximity matrix
  hca <- agnes(diss.matrix, diss = TRUE, method = method)

  
  #cluster membership
  data_clu <- data
  for (i in clu.low:clu.high) {
    clusters <- cutree(hca, i)
    data_clu <-data.frame(data_clu, clusters)
    names(data_clu)[num_var - clu.low + i + 1] <- paste("clu_", i, sep = "" )
  }
  clusters <- data_clu[,(num_var+1):ncol(data_clu)]
  
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
    results <- EVAL(M, clusters, diss = diss.matrix)
    #results <- EVAL(M)
    results1 <- results[[1]]
    results2 <- results[[2]]
  }
  
  clu_results <-  as.list(clusters)
  
  # object for dendrogram creation
  call <- match.call()
  
  dend <- hca[-c(5,6)]
  
  
  if (eval == 1 & prox == 1) {
    list <- list(mem = clu_results, eval = results1, opt = results2, dend = dend, prox = as.dist(diss.matrix), call = call)
  }
  if (eval == 0 & prox == 1) {
    list <- list(mem = clu_results, dend = dend, prox = as.dist(diss.matrix), call = call)
  }
  if (eval == 1 & prox == 0) {
    list <- list(mem = clu_results, eval = results1, opt = results2, dend = dend, call = call)
  }
  if (eval == 0 & prox == 0) {
    list <- list(mem = clu_results, dend = dend, call = call)
  }
  
  attr(list,"class")="nomclust"
  
  return(list)
}

