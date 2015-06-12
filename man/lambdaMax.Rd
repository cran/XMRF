\name{lambdaMax}
\alias{lambdaMax}
\title{Maximum lambda}
\description{Compute the maximum \code{lambda}}
\usage{
lambdaMax(X)
}
\arguments{
  \item{X}{a \emph{n}x\emph{p} data matrix.}
}
\details{Compute the largest value for regularization (maximum \code{lambda}) that gives the null model. 
The maximum \code{lambda} is computed based on the input data matrix, and is the maximum element from column-wise multiplication of data matrix normalized by the
number of observations.
}
\value{a numeric value}
