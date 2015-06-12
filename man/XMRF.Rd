\name{XMRF}
\alias{XMRF}
\title{Markov Random Fields for Exponential Family Distributions}
\description{Infer networks from genomics data using Markov Random Fields specified by node-conditional univariate exponential family distributions.}
\usage{
XMRF(X, method = "LPGM", stability = "bootstrap", N = 100, beta = 0.01, lmin = 0.01, 
   nlams = 20, lambda.path = NULL, parallel = TRUE, nCpus = 4, sym = TRUE, th = 0.01, 
   sth = 0.95, R = max(X), R0 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{X}{a \emph{p}x\emph{n} data matrix.}
	\item{method}{specification of the type of MRF model, default to "\code{LPGM}" for log-linear Poisson-based graphical model. 
		Other allowed methods are "\code{PGM}" for regular Poisson, "\code{TPGM}" for truncated Poisson, "\code{SPGM}" for sublinear Poisson, "\code{GGM}" for Gaussian graphical models, 
		and "\code{ISM}" for Ising model.
	}
  \item{stability}{specification of the stability method, default to "\code{bootstrap}". Another accepted value is "\code{star}" for Stability Approach to Regularization Selection (StARS).}
  \item{N}{number of iterations for stability selection, default to 100.}
  \item{beta}{threshold value on sparsity of the network, default to 0.01.}
  \item{lmin}{ratio of minimum lambda value from the maximum lambda value, default to 0.01.}
  \item{nlams}{number of lambda for regularization, default to 20.}
  \item{lambda.path}{vector lambda used for regularization, default ot \code{NULL}.}
  \item{parallel}{logical value to indicate if the process should be run parallelly in multiple threads, default to \code{TRUE}.}
  \item{nCpus}{number of (maximum) cores to use for parallel execution, default to 4.}
  \item{sym}{logical value to indicate if symmetry is enforced on the inferred edges, default to \code{TRUE}.}
  \item{th}{threshold value for the estimated edge coefficient, default to 0.005.}
  \item{sth}{an inferred edge is retained only if its stability score is greater than \code{sth}, default to 0.9.}
  \item{R}{truncation level for classes "\code{TPGM}" and "\code{SPGM}". The value has to be positive. Default to the maximum value of the input data matrix.}
  \item{R0}{lower-bound truncation level for "\code{SPGM}", default to 0.}
}
\details{
This is the main function of the package that fits exponential family Markov Networks to genomics data. To estimate the network structures using native distribution of the genomics data,
specify the MRF family types in the "\code{method}" parameter. For genomic networks based on next-generation sequencing data, we recommend using the \code{LPGM} family. 
The table at the beginning of the document lists the family type recommended for each of the genomic data platforms.
}
\value{
An object of class \code{GMS} will be returned representing the inferred Markov networks over the regularization path. See \code{\link{GMS}} for details.
}
\references{

Allen, G.I., and Liu, Z. (2012). A Log-Linear graphical model for inferring genetic networks from high-throughput sequencing data. \emph{The IEEE International Conference on Bioinformatics and Biomedicine (BIBM 2012)}.

Allen, G. I., and Liu, Z. (2013). A Local Poisson Graphical Model for Inferring Genetic Networks from Next Generation Sequencing Data. \emph{IEEE Transactions on NanoBioscience}, \bold{12}(3), pp.1-10

Liu, H., Roeder, K., and Wasserman, L. (2010). Stability approach to regularization selection (stars) for high dimensional graphical models. \emph{NIPS 23}, pp.1432?1440.

Meinshausen, N. and Buhlmann, P. (2006). High-dimensional graphs and variable selection with the lasso. \emph{The Annals of Statistics}, \bold{34}(3), pp.1436?1462.

Meinshausen, N. and Buhlmann, P. (2010). Stability selection. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{72}(4), pp.417?473.

Ravikumar, P., Wainwright, M., and Lafferty, J. (2010). High-dimensional ising model selection using l1-regularized logistic regression. \emph{The Annals of Statistics}, \bold{38}(3), pp.1287?1319.

Yang, E., Ravikumar, P.K., Allen, G.I., and Liu, Z. (2012). Graphical models via generalized linear models. \emph{NIPS}, \bold{25}, pp.1367--1375.

Yang, E., Ravikumar, P.K., Allen, G.I., and Liu, Z. (2013a). On graphical models via univariate exponential family distributions. \emph{arXiv preprint arXiv:1301.4183}.

Yang, E., Ravikumar, P.K., Allen, G.I., and Liu, Z. (2013b). On Poisson graphical models. \emph{NIPS}, pp.1718-1726.


}

\seealso{
	\code{\link{XMRF-package}}, \code{\link{GMS}}, \code{\link{plot.GMS}}
}
\examples{
# Example for LPGM
# Refer to the package's introduction for identical example
\dontrun{n = 100}
\dontrun{p = 20}
\dontrun{sim <- XMRF.Sim(n=n, p=p, model="LPGM", graph.type="scale-free")}
\dontrun{simDat <- sim$X}
\dontrun{# Compute the optimal lambda}
\dontrun{lmax = lambdaMax(t(simDat))}
\dontrun{lambda = 0.01* sqrt(log(p)/n) * lmax}
\dontrun{# Run LPGM}
\dontrun{lpgm.fit <- XMRF(simDat, method="LPGM", N=10, lambda.path=lambda)}
\dontrun{ml = plotNet(sim$B, fn="simDat.netPlot.pdf")}
\dontrun{ml = plot(lpgm.fit, fn="lpgm.netPlot_1.pdf", i=1, mylayout=ml)}
\dontrun{plot(lpgm.fit, fn="lpgm.fit.net.pdf")}

}
