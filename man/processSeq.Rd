\name{processSeq}
\alias{processSeq}
\title{Process Sequencing Data for Poisson-based MRFs}
\description{Process and normalize RNA-Sequencing count data into a distribution appropriate for Poisson MRFs.}
\usage{
processSeq(X, quanNorm = 0.75, nLowCount = 20, percentLowCount = 0.95, NumGenes = 500, 
PercentGenes = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{\emph{n}x\emph{p} data matrix.}
  \item{quanNorm}{an optional parameter controlling the quantile for sample normalization, default to 0.75.}
  \item{nLowCount}{minimum read count to decide if to filter a gene, default to 20.}
  \item{percentLowCount}{filter out a gene if it has this percentage of samples less than \code{nLowCount}, default to 0.95.}
  \item{NumGenes}{number of genes to retain in the final data set, default to 500.}
  \item{PercentGenes}{percentage of genes to retain, default to 0.1.}
}
\details{
To process the next-generation sequencing count data into proper distribution (with dispersion removed), the following steps are taken in this function:
	\enumerate{
	\item{Quantile normalization for the samples.}
	\item{Filter out genes with all low counts.}
	\item{Filter genes by maximal variance (if specified).}
	\item{Transform the data to be closer to the Poisson distribution.  A log or power transform is considered and selected based upon the Kolmogorov-Smirnov goodness of fit test.}
	}
}
\value{
a \emph{n} x \code{NumGenes} or \code{PercentGenes} processed data matrix.
}
\examples{
library(XMRF)
data('brcadat')
brca = t(processSeq(t(brcadat), PercentGenes=1))
}
