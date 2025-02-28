% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Subclustering.R
\name{find_subclusters}
\alias{find_subclusters}
\title{Find Sub-clusters in Reference Data}
\usage{
find_subclusters(
  train_object,
  resolution,
  algorithm,
  min_cells,
  verbose = TRUE
)
}
\arguments{
\item{train_object}{Training object}

\item{resolution}{Clustering resolution}

\item{algorithm}{Clustering algorithm}

\item{min_cells}{Minimum cells per cluster}

\item{verbose}{Show progress}
}
\value{
Vector of sub-cluster assignments
}
\description{
Find Sub-clusters in Reference Data
}
\keyword{internal}
