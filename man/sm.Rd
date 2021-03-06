% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sm.R
\name{sm}
\alias{sm}
\title{Simple Matching Coefficient (SM)}
\usage{
sm(data)
}
\arguments{
\item{data}{A data.frame or a matrix with cases in rows and variables in colums.}
}
\value{
The function returns an object of the class "dist".
\cr
}
\description{
The function calculates a dissimilarity matrix based on the SM similarity measure.
}
\details{
The simple matching coefficient (Sokal, 1958) represents the simplest way of measuring similarity. It does not impose any weights.
By a given variable, it assigns the value 1 in case of match and value 0 otherwise.
}
\examples{
# sample data
data(data20)

# dissimilarity matrix calculation
prox.sm <- sm(data20)

}
\references{
Boriah S., Chandola V., Kumar V. (2008). Similarity measures for categorical data: A comparative evaluation.
In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254.
 \cr
 \cr
Sokal R., Michener C. (1958). A statistical method for evaluating systematic relationships. In: Science bulletin, 38(22),
The University of Kansas.  
 \cr
}
\seealso{
\code{\link[nomclust]{eskin}},
\code{\link[nomclust]{good1}},
\code{\link[nomclust]{good2}},
\code{\link[nomclust]{good3}},
\code{\link[nomclust]{good4}},
\code{\link[nomclust]{iof}},
\code{\link[nomclust]{lin}},
\code{\link[nomclust]{lin1}},
\code{\link[nomclust]{of}},
\code{\link[nomclust]{ve}},
\code{\link[nomclust]{vm}}.
}
\author{
Zdenek Sulc. \cr Contact: \email{zdenek.sulc@vse.cz}
}
