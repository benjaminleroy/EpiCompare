% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bands.R
\name{check_inside_elipsoid}
\alias{check_inside_elipsoid}
\title{assert if observation is inside elipsoid}
\usage{
check_inside_elipsoid(data, Sigma, center, bound, suppress_warning = FALSE)
}
\arguments{
\item{data}{data.frame or matrix of data (row is observation), ncol = p,
scalar columns}

\item{Sigma}{covariance matrix (p x p)}

\item{center}{center of elipsoid (vector length p)}

\item{bound}{contraint for equation of ellipsoid}

\item{suppress_warning}{logic to suppress warning if just returning all
\code{FALSE} relative to non PSD Sigma or bound <= 0.}
}
\value{
boolean vector if data is contained in ellipsoid. All \code{FALSE} if
  Sigma not PSD
}
\description{
See
https://stats.stackexchange.com/questions/29860/confidence-interval-of-multivariate-gaussian-distribution
for details.
}
