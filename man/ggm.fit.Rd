\name{ggm.fit}
\alias{ggm.fit}
\docType{data}
\title{Fitted Gaussian Graphical Models}
\description{
An example fitted Gaussian graphical model 
}
\usage{data("ggm.fit")}
\format{
A \code{\link{GMS}} object.
}
\details{
This fitted Gaussian graphical model is included as the example model for demonstrating the usage of other functions. 
The model was fitted with a simulated multivariate Gaussian data of scale-free graph structure, 30 nodes, and 200 observations. 
StARS stability selection of 100 iterations was used to identify the optimal network over the regularization path of 10 parameters. 
}

\keyword{datasets}
\seealso{
	\code{\link{GMS}}, \code{\link{plot.GMS}}
}