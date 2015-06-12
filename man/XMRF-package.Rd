\name{XMRF-package}
\alias{XMRF-package}
\alias{XMRF-package}
\docType{package}
\title{
A R Package to Fit Markov Networks to High-throughput Genomics Data
}
\description{
An R package to learn and visualize the underlying relationships between genes from various types of high-throughput genomics data.
}
\details{
\tabular{rl}{
Package: \tab XMRF\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-06-12\cr
License: \tab GPL-2\cr
}


Technological advances have produced large amounts of high-throughput "omics" data that allow us to study the complicated interactions between genes, mutations, aberrations, and epi-genetic markers. 
Markov Random Fields (MRFs), or Markov Networks, enable us to estimate these genomics networks via sparse, high-dimensional undirected graphical models. 

\bold{}

Here, we provide the community a convenient and useful tool to learn the complex genomics networks from various types of high-throughput genomics data. 
This package encodes the recently proposed parametric family of graphical models based on node-conditional univariate exponential family distributions (Yang \emph{et. al}, 2012, 2013a). 
Specifically, our package has methods for estimating Gaussian graphical model (Meinshausen and Buhlmann, 2006), Ising model (Ravikumar \emph{et. al}, 2010), and the Poisson family graphical models (Allen and Liu, 2012, 2013; Yang \emph{et. al} 2013b).
These models can be used to estimate genetic networks from a variety of data types:

\tabular{lll}{
		Genomics Data \tab Type \tab XMRF Family \cr 
		======================== \tab ========== \tab ============\cr
		RNA-Seq or miRNA-Seq \tab Counts \tab \code{LPGM}, \code{SPGM} \cr
		Microarray or Methylation array \tab Continuous \tab \code{GGM} \cr
		Mutation or CNVs \tab Binary \tab \code{ISM} \cr
}


To estimate the network structures from different types of genomics data with this package, users simply need to specify the \code{"method"} in the main function, 
for example \code{XMRF(..., method="LPGM")} to fit LPGM to next-generation sequencing data.

\bold{}

In this package, we implement the neighborhood selection graph estimation technique by proximal or projected gradient descent using warm starts over the range of regularization parameters. 
This technique allows estimation of the neighborhood for each node independently and thus can be done in parallel, thus speeding computation.

\bold{}

This package also implements two data-driven approaches to select the sparsity of a fitted network: a stability-based approach for a single regularization value over many bootstrap resamples (Meinshausen and Buhlmann, 2010), or
computed over a range of regularization values with StARS ({Liu \emph{et. al.}, 2010}).

}
\author{
Ying-Wooi Wan, Genevara Allen, Yulia Baker, Eunho Yang, Pradeep Ravikumar, and Zhandong Liu

Maintainer: Ying-Wooi Wan<yingwoow@bcm.edu>
}
\references{
Allen, G.I., and Liu, Z. (2012). A Log-Linear graphical model for inferring genetic networks from high-throughput sequencing data. \emph{The IEEE International Conference on Bioinformatics and Biomedicine (BIBM 2012)}.

Allen, G. I., and Liu, Z. (2013). A Local Poisson Graphical Model for Inferring Genetic Networks from Next Generation Sequencing Data. \emph{IEEE Transactions on NanoBioscience}, \bold{12}(3), pp.1-10

Liu, H., Roeder, K., and Wasserman, L. (2010). Stability approach to regularization selection (stars) for high dimensional graphical models. \emph{NIPS 23}, pp.1432-1440.

Meinshausen, N. and Buhlmann, P. (2006). High-dimensional graphs and variable selection with the lasso. \emph{The Annals of Statistics}, \bold{34}(3), pp.1436-1462.

Meinshausen, N. and Buhlmann, P. (2010). Stability selection. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{72}(4), pp.417-473.

Ravikumar, P., Wainwright, M., and Lafferty, J. (2010). High-dimensional ising model selection using l1-regularized logistic regression. \emph{The Annals of Statistics}, \bold{38}(3), pp.1287-1319.

Yang, E., Ravikumar, P.K., Allen, G.I., and Liu, Z. (2012). Graphical models via generalized linear models. \emph{NIPS}, \bold{25}, pp.1367--1375.

Yang, E., Ravikumar, P.K., Allen, G.I., and Liu, Z. (2013a). On graphical models via univariate exponential family distributions. \emph{arXiv preprint arXiv:1301.4183}.

Yang, E., Ravikumar, P.K., Allen, G.I., and Liu, Z. (2013b). On Poisson graphical models. \emph{NIPS}, pp.1718-1726.

}
\keyword{ package }
\seealso{
	\code{\link{XMRF}}
}
\examples{
	library(XMRF)
	
	n = 100
	p = 20
	sim <- XMRF.Sim(n=n, p=p, model="LPGM", graph.type="scale-free")
	simDat <- sim$X
	
	# Compute the optimal lambda
	lmax = lambdaMax(t(simDat))
	lambda = 0.01* sqrt(log(p)/n) * lmax
	# Run LPGM
	lpgm.fit <- XMRF(simDat, method="LPGM", N=10, lambda.path=lambda, parallel=FALSE)
	
	# Print the fitted Markov networks
	lpgm.fit
	
	ml = plotNet(sim$B)
	ml = plot(lpgm.fit, mylayout=ml)
	}
