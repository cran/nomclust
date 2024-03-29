% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nomclust.R
\name{nomclust}
\alias{nomclust}
\title{Hierarchical Clustering of Nominal Data}
\usage{
nomclust(
  data,
  measure = "lin",
  method = "average",
  clu.high = 6,
  eval = TRUE,
  prox = 100,
  var.weights = NULL
)
}
\arguments{
\item{data}{A data.frame or a matrix with cases in rows and variables in columns.}

\item{measure}{A character string defining the similarity measure used for computation of proximity matrix in HCA:
\code{"anderberg"}, \code{"burnaby"}, \code{"eskin"}, \code{"gambaryan"}, \code{"goodall1"}, \code{"goodall2"}, \code{"goodall3"}, \code{"goodall4"}, \code{"iof"}, \code{"lin"}, \code{"lin1"}, \code{"of"}, \code{"sm"}, \code{"smirnov"}, \code{"ve"}, \code{"vm"}.}

\item{method}{A character string defining the clustering method. The following methods can be used: \code{"average"}, \code{"complete"}, \code{"single"}.}

\item{clu.high}{A numeric value expressing the maximal number of cluster for which the cluster memberships variables are produced.}

\item{eval}{A logical operator; if TRUE, evaluation of the clustering results is performed.}

\item{prox}{A logical operator or a numeric value. If a logical value TRUE indicates that the proximity matrix is a part of the output. A numeric value (integer) of this argument indicates the maximal number of cases in a dataset for which a proximity matrix will occur in the output.}

\item{var.weights}{A numeric vector setting weights to the used variables. One can choose the real numbers from zero to one.}
}
\value{
The function returns a list with up to six components.
\cr
\cr
The \code{mem} component contains cluster membership partitions for the selected numbers of clusters in the form of a list.
\cr
\cr
The \code{eval} component contains up to 13 evaluation criteria as vectors in a list. Namely, Within-cluster mutability coefficient (WCM), Within-cluster entropy coefficient (WCE),
Pseudo F Indices based on the mutability (PSFM) and the entropy (PSFE), Bayesian (BIC), and Akaike (AIC) information criteria for categorical data, the BK index, Category Utility (CU), Category Information (CI), Hartigan Mutability (HM), Hartigan Entropy (HE) and, if the prox component is present, the silhouette index (SI) and the Dunn index (DI).
\cr
\cr
The \code{opt} component is present in the output together with the \code{eval} component. It displays the optimal number of clusters for the evaluation criteria from the \code{eval} component, except for WCM and WCE, where the optimal number of clusters is based on the elbow method.
\cr
\cr
The \code{dend} component can be found in the output together with the \code{prox} component. It contains all the necessary information for dendrogram creation.
\cr
\cr
The \code{prox} component contains the dissimilarity matrix in the form of the "dist" object.
\cr
\cr
The \code{call} component contains the function call.
}
\description{
The function performs and evaluates hierarchical cluster analysis of nominal data.
}
\details{
The function runs hierarchical cluster analysis (HCA) with objects characterized by nominal variables (without natural order of categories).
 It completely covers the clustering process, from the dissimilarity matrix calculation to the cluster quality evaluation. The function enables a user to choose from the similarity measures for nominal data summarized by (Boriah et al., 2008) and by (Sulc and Rezankova, 2019). 
 Next, it offers to choose from three linkage methods that can be used for categorical data. It is also possible to assign user-defined variable weights. The obtained clusters can be evaluated by up to 13 evaluation criteria (Sulc et al., 2018) and (Corter and Gluck, 1992). The output of the nomclust() function may serve as an input for the visualization functions \emph{dend.plot} and \emph{eval.plot} in the nomclust package.
}
\examples{
# sample data
data(data20)

# creating an object with results of hierarchical clustering of 
hca.object <- nomclust(data20, measure = "lin", method = "average",
 clu.high = 5, prox = TRUE)

# assigning variable weights
hca.weights <- nomclust(data20, measure = "lin", method = "average",
 clu.high = 5, prox = TRUE, var.weights = c(0.7, 1, 0.9, 0.5, 0))

# quick clustering summary
summary(hca.object)

# quick cluster quality evaluation
print(hca.object)

# visualization of the evaluation criteria
eval.plot(hca.object)

# a quick dendrogram
plot(hca.object)

# a dendrogram with three designated clusters
dend.plot(hca.object, clusters = 3)

# obtaining values of evaluation indices as a data.frame
data20.eval <- as.data.frame(hca.object$eval)

# getting the optimal numbers of clusters as a data.frame
data20.opt <- as.data.frame(hca.object$opt)

# extracting cluster membership variables as a data.frame
data20.mem <- as.data.frame(hca.object$mem)

# obtaining a proximity matrix
data20.prox <- as.matrix(hca.object$prox)

# setting the maximal number of objects for which a proximity matrix is provided in the output to 30
hca.object <- nomclust(data20, measure = "iof", method = "complete",
 clu.high = 5, prox = 30)
 
# transforming the nomclust object to the class "hclust"
hca.object.hclust <- as.hclust(hca.object)

# transforming the nomclust object to the class "agnes, twins"
hca.object.agnes <- as.agnes(hca.object)


}
\references{
Boriah S., Chandola V. and Kumar, V. (2008). Similarity measures for categorical data: A comparative evaluation.
In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
\cr
\cr
Corter J.E., Gluck M.A. (1992). Explaining basic categories: Feature predictability and information. Psychological Bulletin 111(2), p. 291–303.
\cr
\cr
Sulc Z., Cibulkova J., Prochazka J., Rezankova H. (2018). Internal Evaluation Criteria for Categorical Data in Hierarchical Clustering: Optimal Number of Clusters Determination, Metodoloski Zveski, 15(2), p. 1-20.
\cr
\cr
Sulc Z. and Rezankova H. (2019). Comparison of Similarity Measures for Categorical Data in Hierarchical Clustering. Journal of Classification, 35(1), p. 58-72. DOI: 10.1007/s00357-019-09317-5.
}
\seealso{
\code{\link[nomclust]{evalclust}}, \code{\link[nomclust]{nomprox}}, \code{\link[nomclust]{eval.plot}}, \code{\link[nomclust]{dend.plot}}.
}
\author{
Zdenek Sulc. \cr Contact: \email{zdenek.sulc@vse.cz}
}
