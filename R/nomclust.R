#' Hierarchical Cluster Analysis for Nominal Data
#' 
#' @description The \code{nomclust()} function runs hierarchical cluster analysis (HCA) with objects characterized by nominal (categorical) variables. It completely covers the clustering process, from the proximity matrix calculation to the evaluation of the quality of clustering.
#' The function contains thirteen similarity measures for nominal data summarized in (Boriah et al., 2008) or introduced by Morlini and Zani in (Morlini and Zani, 2012), and by (Sulc and Rezankova, 2019). 
#' It offers three linkage methods that can be used for categorical data. The obtained clusters can be evaluated by seven evaluation criteria, see (Sulc et al., 2018). The output of the \code{nomclust()} function may serve as an input for visualization functions in the \bold{nomclust} package.
#' 
#' 
#' @param data A \emph{data.frame} or a \emph{matrix} with cases in rows and variables in colums.
#' 
#' @param measure A \emph{character} string defining the similarity measure used for computation of proximity matrix in HCA:
#' \code{"eskin"}, \code{"good1"}, \code{"good2"}, \code{"good3"}, \code{"good4"}, \code{"iof"}, \code{"lin"}, \code{"lin1"}, \code{"morlini"}, \code{"of"}, \code{"sm"}, \code{"ve"}, \code{"vm"}.
#' 
#' @param method A \emph{character} string defining the clustering method. The following methods can be used: \code{"average"}, \code{"complete"}, \code{"single"}.
#' 
#' @param clu.high  A \emph{numeric} value expressing the maximal number of cluster for which the cluster memberships variables are produced.
#' 
#' @param eval A \emph{logical} operator; if TRUE, evaluation of the clustering results is performed.
#' 
#' @param prox A \emph{logical} operator or a numeric value. If a logical value TRUE indicates that the proximity matrix is a part of the output. A numeric value (integer) of this argument indicates the maximal number of cases in a dataset for which a proximity matrix will occur in the output.
#' 
#' @param opt A \emph{logical} operator; if TRUE, the time optimization method is run to substantially decrease computation time of the dissimilarity matrix calcation. Time optimalization method cannot be run if the proximity matrix is to be produced. In such a case, this parameter is automatically set to FALSE.
#' 
#' @return The function returns a \emph{list} with up to five components.
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
#' \cr
#' \cr
#' The \code{prox} component contains the dissimilarity matrix in a form of a \emph{matrix}.
#' \cr
#' \cr
#' The \code{dend} component can be found in the output only together with the \code{prox} component. It contains all the necessary information for dendrogram creation.
#' 
#'
#'@references
#' Boriah S., Chandola V. and Kumar, V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
#' \cr
#' \cr
#' Morlini I. and Zani S. (2012). A new class of weighted similarity indices using polytomous variables. Journal of Classification, 29(2), p. 199-226.
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
#'  clu.high = 5, prox = TRUE, opt = FALSE)
#' 
#' # obtaining values of evaluation indices
#' data20.eval <- hca.object$eval
#' 
#' # getting the optimal numbers of clusters
#' data20.opt <- hca.object$opt
#' 
#' # extracting cluster membership variables
#' data20.mem <- hca.object$mem
#' 
#' # extracting cluster membership variables as a data frame
#' data20.mem <- as.data.frame(hca.object$mem)
#' 
#' # obtaining a proximity matrix
#' data20.prox <- hca.object$prox
#' 
#' # setting the maximal number of objects for which a proximity matrix is provided in the output to 30
#' hca.object <- nomclust(data20, measure = "lin", method = "average",
#'  clu.high = 5, prox = 30, opt = FALSE)
#' 
#' # generating of a larger dataset containing repeatedly occuring objects
#' set.seed(150)
#' sample150 <- sample(1:nrow(data20), 150, replace = TRUE)
#' data150 <- data20[sample150, ]
#' 
#' # running hierarchical clustering WITH the time optimization
#' start <- Sys.time()
#' hca.object.opt.T <- nomclust(data150, measure = "lin", opt = TRUE)
#' end <- Sys.time()
#' end - start
#' 
#' # running hierarchical clustering WITHOUT the time optimization
#' start <- Sys.time()
#' hca.object.opt.F <- nomclust(data150, measure = "lin", opt = FALSE)
#' end <- Sys.time()
#' end - start
#' 
#' @export

nomclust <- function (data, measure = "lin", method = "average", clu.high = 6, eval = TRUE, 
                       prox = 100, opt = TRUE) {
  
  clu.low = 2
  
  # elimination of other clustering methods
  if (method %in% c("single", "average", "complete") == FALSE) {
    stop("An invalid clustering method was chosen.")
  }
  
  # number of clusters cannot exceed the parameter clu.high
  if (nrow(data)<clu.high) {
    stop("The 'clu.high' argument cannot exceed the number of clustered objects.")
  }
  
  # calculate proximity matrix for up to maximal size of a dataset
  if (prox != FALSE & prox != TRUE) {
    if (is.numeric(prox) == T) {
      if (nrow(data) <= abs(as.integer(prox))) {
        prox <- T
        opt <- F
      } else
        prox <- F
    } else
      stop("The 'prox' argument should be of the 'numeric' type.")
  }
  

  # to prevent publication of a proximity matrix with unordered rows
  if (opt == TRUE & prox == TRUE) {
    opt <- FALSE
    print("The time optimization method cannot be run if the proximity matrix is to be produced. The standard calculation method is used instead.")
  }
  
  # for which measures the optimizalization will be used
  if (opt == TRUE) {
    if (measure %in% c("eskin", "good1", "good2", "good3", "good4", "iof", "of", "lin", "lin1", "sm", "ve", "vm")) {
      opt <- TRUE
      # print("Time optimalization in process.")
    } else if (measure %in% c("morlini")) {
      opt <- FALSE
      print("The time optimization cannot be used with the 'morlini' similarity measure. The standard calculation method is used instead.")
    } else {
      stop("Invalid name of the similarity measure.")
    }
  } else {
      if (measure %in% c("eskin", "good1", "good2", "good3", "good4", "iof", "of", "lin", "lin1", "morlini", "sm", "ve", "vm")) {
        opt <- FALSE
        # print("Standard calculation method in process.")
      } else {
        stop("An invalid name of the similarity measure was used.")
      }
  }
  
  
  measure <- paste(measure, "_fx", sep = "")
  
  sim_measure <- get(measure)
  
  rnames <- row.names(data)
  
  # if matrix, coerce to data.frame
  if(is.matrix(data) == 1) {
    data <- as.data.frame(data)
  }
  
  # recoding everything to factors and then to numeric values
  indx <- sapply(data, is.factor)
  data[!indx] <- sapply(data[!indx], function(x) as.factor(x))
  data <- as.data.frame(unclass(data))
  data <- sapply(data, function(x) as.numeric(x))
  data <- as.data.frame(data)
  
  
  # number of variables of dataset
  num_var <- ncol(data)
  
  # max number of categories
  num_cat <- sapply(data, function(x) length(unique(x)))
  max_num_cat <- max(num_cat)
  
  # frequencies of categories in all variables
  freq.table <- freq.abs(data)
  
  
  if (opt == TRUE) {
    
    # pridat novy sloupec s poradim hodnot v pripade, ze rownames jsou slova
    # ordering the data according to all variables
    
    row.names(data) <- NULL
    
    for (i in ncol(data):1 ) {       
      data <- data[order(data[,i]), ]
    }
    
    # adding vector keeping the original order of cases
    order <- as.numeric(row.names(data))
    
    # library(plyr)
    # data reduction
    data_rw <- ddply(data, names(data), as.data.frame(nrow))
    data_r <- data_rw[,-ncol(data_rw)] # reduced
    data_w <- data_rw[,ncol(data_rw)]  # weight
    
    # VYRESIT - ASI PODMINKOU
    # if (measure %in% c("morlini")) {
    #   prox_matrix <- sim_measure(data = data, freq = freq.table) # for morlini
    # } else {
      prox_matrix <- sim_measure(data = data_r, freq = freq.table)
    # }
    
    
    prox2 <-  matrix(data=0,nrow=nrow(data),ncol=nrow(data))
    
    index1 <- 1
    index2 <- data_w[1]
    
    for (i in 1:nrow(prox_matrix)) {
      index3 <- 1
      index4 <- data_w[1]
      for (j in 1:ncol(prox_matrix)) {
        prox2[index1:index2,index3:index4] <- prox_matrix[i,j]
        index3 <- index4 + 1
        index4 <- index4 + data_w[j+1]
      }
      index1 <- index2 + 1
      index2 <- index2 + data_w[i+1]
    }
    
    prox_matrix <- prox2
    
    
  } else {
    #computing the proximity matrices
    prox_matrix <- sim_measure(data, freq = freq.table)
    row.names(prox_matrix) <- rnames
  }
  
  #row.names(prox_matrix) <- rnames VYPNUTE DO DOBY NEZ BUDE OPRAVENE PORADI HODNOT V PRIPADE OPTIMALIZACE
  
  
  
  #hierarchical cluster analysis, where "prox" is a proximity matrix
  hca <- agnes(prox_matrix, diss = TRUE, method = method)

  #cluster membership
  data_clu <- vector(mode="numeric", length = nrow(data))
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
  
  # obtaining the original cluster order
  if (opt == TRUE) {
    clusters <- clusters[order(as.numeric(row.names(clusters))), ]
  }
  

  clu_results <-  as.list(clusters)
  
  
  # object for dendrogram creation
  dend <- hca[-c(5,6)]
  
  # obtaining a correct labels in a dendrogram
  # if (opt == T) {
  #   order.labels <- cbind(order, dend$order.lab)
  #   order.labels <- order.labels[order(as.numeric(order.labels[,1])), ]
  #   dend$order.lab <- order.labels[,2]
  # }
  
  if (eval == 1 & prox == 1 & opt == 0) {
    list <- list(mem = clu_results, eval = results1, opt = results2, prox = prox_matrix, dend = dend)
  }
  if (eval == 0 & prox == 1 & opt == 0) {
    list <- list(mem = clu_results, prox = prox_matrix, dend = dend)
  }
  if (eval == 1 & prox == 0 & opt == 0) {
    list <- list(mem = clu_results, eval = results1, opt = results2, dend = dend)
  }
  if (eval == 0 & prox == 0 & opt == 0) {
    list <- list(mem = clu_results, dend = dend)
  }
  if (eval == 1 & prox == 0 & opt == 1) {
    list <- list(mem = clu_results, eval = results1, opt = results2)
  }
  if (eval == 0 & prox == 0 & opt == 1) {
    list <- list(mem = clu_results)
  }
  
  return(list)
}
