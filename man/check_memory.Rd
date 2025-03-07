% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_memory}
\alias{check_memory}
\title{Memory Usage Reporter}
\usage{
check_memory(threshold = 10000)
}
\arguments{
\item{threshold}{Memory threshold in MB to trigger warning}
}
\value{
Named vector of memory statistics
}
\description{
Reports current memory usage with garbage collection
}
\keyword{internal}
