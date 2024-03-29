% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evalclust.R
\name{evalclust}
\alias{evalclust}
\title{Cluster Quality Evaluation of Nominal Data Hierarchical Clustering}
\usage{
evalclust(data, clusters, diss = NULL)
}
\arguments{
\item{data}{A data.frame or a matrix with cases in rows and variables in columns.}

\item{clusters}{A data.frame or a list of cluster memberships obtained based on the dataset defined in the parameter \code{data} in the form of a sequence from the two-cluster solution to the maximal-cluster solution.}

\item{diss}{An optional parameter. A matrix or a dist object containing dissimilarities calculated based on the dataset defined in the parameter \code{data}.}
}
\value{
The function returns a list with three components.
\cr
\cr
The \code{eval} component contains up to 13 evaluation criteria as vectors in a list. Namely, Within-cluster mutability coefficient (WCM), Within-cluster entropy coefficient (WCE),
Pseudo F Indices based on the mutability (PSFM) and the entropy (PSFE), Bayesian (BIC), and Akaike (AIC) information criteria for categorical data, the BK index, Category Utility (CU), Category Information (CI), Hartigan Mutability (HM), Hartigan Entropy (HE) and, if the prox component is present, the silhouette index (SI) and the Dunn index (DI).
\cr
\cr
The \code{opt} component is present in the output together with the \code{eval} component. It displays the optimal number of clusters for the evaluation criteria from the \code{eval} component, except for WCM and WCE, where the optimal number of clusters is based on the elbow method.
\cr
\cr
The \code{call} component contains the function call.
}
\description{
The function evaluates clustering results by a set of evaluation criteria (cluster validity indices).
}
\details{
The function calculates a set of evaluation criteria if the original dataset and the cluster membership variables are provided. 
 The function calculates up to 13 evaluation criteria described by (Sulc et al., 2018) and (Corter and Gluck, 1992) and provides the optimal number of clusters based on these criteria. 
 It is primarily focused on evaluating hierarchical clustering results obtained by similarity measures different from those that occur in the nomclust package. 
 Thus, it can serve for the comparison of various similarity measures for categorical data.
}
\examples{
# sample data
data(data20)

# creating an object with results of hierarchical clustering
hca.object <- nomclust(data20, measure = "iof", method = "average", clu.high = 7)

# the cluster memberships
data20.clu <- hca.object$mem

# obtaining evaluation criteria for the provided dataset and cluster memberships
data20.eval <- evalclust(data20, clusters = data20.clu)

# visualization of the evaluation criteria
eval.plot(data20.eval)

# silhouette index can be calculated if the dissimilarity matrix is provided
data20.eval <- evalclust(data20, clusters = data20.clu, diss = hca.object$prox)
eval.plot(data20.eval, criteria = "SI")

}
\references{
Corter J.E., Gluck M.A. (1992). Explaining basic categories: Feature predictability and information. Psychological Bulletin 111(2), p. 291–303.
\cr
\cr
Sulc Z., Cibulkova J., Prochazka J., Rezankova H. (2018). Internal Evaluation Criteria for Categorical Data in Hierarchical Clustering: Optimal Number of Clusters Determination, Metodoloski Zveski, 15(2), p. 1-20.
}
\seealso{
\code{\link[nomclust]{nomclust}}, \code{\link[nomclust]{nomprox}}, \code{\link[nomclust]{eval.plot}}.
}
\author{
Zdenek Sulc. \cr Contact: \email{zdenek.sulc@vse.cz}
}
