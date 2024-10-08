% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_type.R
\name{get_relevant_ids}
\alias{get_relevant_ids}
\title{Function to get relevant identifiers from the environment}
\usage{
get_relevant_ids(identifiers, dt)
}
\arguments{
\item{identifiers}{A character vector of identifier names to fetch from the environment}

\item{dt}{A data.table containing the columns to be checked against the identifiers}
}
\value{
A character vector of relevant identifiers that are present in the data.table
}
\description{
Function to get relevant identifiers from the environment
}
\keyword{internal}
