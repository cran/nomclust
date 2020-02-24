\name{CA.methods}
\docType{data}
\alias{CA.methods}
\title{Selected clustering algorithms}
\description{
  The dataset contains five different characteristics of 24 clustering algorithms. The "Type" variable expresses the principle on which the clustering is based. There are five possible categories: density, grid, hierarchical, model-based, and partitioning. The binary variable "OptClu" indicates if the clustering algorithm offers the optimal number of clusters. The variable "Large" indicates if the clustering algorithm was designed to cluster large datasets. The "TypicalType" variable presents the typical data type for which the clustering algorithm was determined. There are three possible categories: categorical, mixed, and quantitative. Since some clustering algorithms support more data types, the binary variable "MoreTypes" indicates this support.
  }
\usage{data("CA.methods")}
\format{A data frame containing 5 variables and 24 cases.}
\source{created by the authors of the nomclust package}
\keyword{clustering}