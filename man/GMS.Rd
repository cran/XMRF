\name{GMS}
\alias{GMS}
\title{Graphical Models Object}
\description{
This class of objects is returned by the \code{XMRF} function in this package, to represent the fitted Markov Network(s).
Objects of this class have the \code{print} method to display the core information of the fitted models and \code{plot} method to plot the optimal Markov Network.
}
\arguments{
  \item{v}{vector of (\code{nlams} length) network variability measured for each regularization level.}
  \item{D}{a list of \emph{p}x\emph{p} matrices of stability scores of inferred edges of each network along the regularization path.}
  \item{lambda.path}{numeric vector used for regularization path.}
  \item{opt.lambda}{lambda value that gives the optimal network (network with maximum variability).}
  \item{network}{a list of \emph{p}x\emph{p} matrices of fitted networks along the regularization path.}
  \item{opt.index}{index of the regularization value that gives the optimal network.}
}
\seealso{
	\code{\link{XMRF}}, \code{\link{plot.GMS}}
}
