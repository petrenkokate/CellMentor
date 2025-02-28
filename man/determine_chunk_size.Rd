% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_csfnmf.R, R/train_functions.R
\name{determine_chunk_size}
\alias{determine_chunk_size}
\title{Determine optimal chunk size}
\usage{
determine_chunk_size(n_cells, available_memory = 1000)

determine_chunk_size(n_cells, available_memory = 1000)
}
\arguments{
\item{n_cells}{Number of cells}

\item{available_memory}{Available memory in MB}
}
\value{
Optimal chunk size

Optimal chunk size
}
\description{
Determine optimal chunk size

Helper function for optimal chunk size determination
}
\keyword{internal}
