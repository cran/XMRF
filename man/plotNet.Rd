\name{plotNet}
\alias{plotNet}
\title{Plot Network}
\description{Plot a network with specific layout.}
\usage{
plotNet(net, fn = "", th = 1e-06, mylayout = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{a square adjacency matrix of the network to be plotted.}
  \item{fn}{file name to save the network plot. Default to be an empty string, 
  so the network is plotted to the standard output (screen). NOTE: if a file name is specified, it should be file name for PDF file.}
  \item{th}{numeric value, default to 1e-06. To specify the threshold if the estimated coefficient between two variables is to be considered connected.}
  \item{mylayout}{graph layout to draw the network, default to \code{NULL}.}
}
\details{This function serves as the alternative plotting function to allow users to plot a specific network with specific layout, such as plotting the simulated network.}
\value{
Returns the layout object from \code{igraph} package - numeric matrix of two columns and the rows with the same number as the number of vertices. }
\examples{
library(XMRF)
n = 200
p = 30
sim <- XMRF.Sim(n=n, p=p, model="LPGM", graph.type="scale-free")
ml = plotNet(sim$B)
}
