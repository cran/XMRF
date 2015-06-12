\name{plot.GMS}
\alias{plot.GMS}
\alias{plot}
\title{Plot GMS Object}
\description{
	Default function to plot the network of a GMS object. 
}
\usage{
\method{plot}{GMS}(x, fn = "", th = 1e-06, i = NULL, mylayout = NULL, vars = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{GMS} object.}
  \item{fn}{file name to save the network plot; default to be an empty string, 
  so the network is plotted to the standard output (screen). NOTE: if a file name is specified, it should be file name for PDF file.}
  \item{th}{numeric value, default to 1e-06. To specify the threshold if the estimated coefficient between two variables is to be considered connected.}
  \item{i}{index of the network (along the regularization path) to be plotted. Default to \code{NULL} to plot the optimal network.}
  \item{mylayout}{graph layout to draw the network, default to \code{NULL}.}
  \item{vars}{vector of variable names, default to \code{NULL}.}
  \item{\dots}{other generic arguments for \code{plot} method.}
}
\details{
	This is the default plotting function for \code{GMS} objects (Markov Networks inferred over a regularization path). Refer to \code{\link{GMS}} for details on \code{GMS} object.
	This function will plot the optimal network on the screen by default.  However, given a file name, the plot will be saved to a PDF file. 
	Also, given a specific index corresponding to the index of \code{lambda.path}, the associated network will be plotted.
	{}
	The network will be plotted in force-directed layout (\code{layout.fruchterman.reingold} with default parameters implemented in \code{igraph} package).
}
\value{
	Returns the layout object from \code{igraph} package - numeric matrix of two columns and the rows with the same number as the number of vertices. 
}
\seealso{
	\code{\link{GMS}}
}
\examples{
library(XMRF)
data('ggm.fit')
plot(ggm.fit, fn="ggm.fit.net.pdf")
}
