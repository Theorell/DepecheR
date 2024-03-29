% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DepecheR.R
\docType{package}
\name{DepecheR-package}
\alias{DepecheR}
\alias{DepecheR-package}
\title{DepecheR: Determination of essential phenotypic elements of clusters in high-dimensional entities}
\description{
The purpose of this package is to identify traits in a dataset that can separate groups. This is done on two levels. First, clustering is performed, using an implementation of sparse K-means. Secondly, the generated clusters are used to predict outcomes of groups of individuals based on their distribution of observations in the different clusters. As certain clusters with separating information will be identified, and these clusters are defined by a sparse number of variables, this method can reduce the complexity of data, to only emphasize the data that actually matters.
}
\details{
The package is indirectly clearly dependent on Rtsne for generation of
output. See comment on the github wiki for how to speed Rtsne up.
}
\author{
\strong{Maintainer}: Jakob Theorell \email{jakob.theorell@ki.se}

Authors:
\itemize{
  \item Axel Theorell \email{axel.theorell@bsse.ethz.ch}
}

}
