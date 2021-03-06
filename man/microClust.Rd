% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microClust.R
\name{microClust}
\alias{microClust}
\title{This function is the core of the neighSmooth. See the documentation
there for details.}
\usage{
microClust(
  dataCenter,
  dataNeigh,
  dataReturn,
  method = "median",
  k = 11,
  trim = 0
)
}
\arguments{
\item{dataCenter}{The original data.}

\item{dataNeigh}{The data for the neighbors. Often stronly overlapping with
the dataCenter, but for internal reasons, this data cloud is larger
than the dataCenter cloud.}

\item{dataReturn}{The neighbor data that should be aggregated and sent back.}

\item{method}{Should median or mean be calculated?}

\item{k}{Number of neighbors.}

\item{trim}{If mean of the neighbors is returned, should it be calculated
with trimming?}
}
\value{
A dataset with the same shape as dataCenter, filled with
aggregated information from the k nearest neighbors.
}
\description{
This function is the core of the neighSmooth. See the documentation
there for details.
}
\keyword{internal}
