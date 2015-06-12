\name{plotGML}
\alias{plotGML}
\title{Plot Network in GML}
\description{Plot the network in graph modeling language (GML).}
\usage{
plotGML(x, fn = "", th = 1e-06, i = NULL, weight = FALSE, vars = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{GMS} object.}
  \item{fn}{file name to save the GML file.}
  \item{th}{numeric value, default to 1e-06. To specify the threshold if the estimated coefficient between two variables is to be considered connected.}
  \item{i}{index of the network (along the regularization path) to be plotted. Default to \code{NULL} for optimal network.}
  \item{weight}{boolean value to indicate if writing the stability on the inferred edges, default to \code{FALSE}.}
  \item{vars}{vector of variable names, default to \code{NULL}.}
}
