\name{XMRF.Sim}
\alias{XMRF.Sim}
\title{Generate simulated data from XMRF models}
\description{Generate data from different multivariate distributions with different network structures.}
\usage{
XMRF.Sim(n = 100, p = 50, model = "LPGM", graph.type = "scale-free")
}
\arguments{
  \item{n}{number of samples, default to 100.}
  \item{p}{number of variables, default to 50.}
  \item{model}{Markov Network models to indicate the distribution family of the data to be generated, default to "\code{LPGM}". 
		Other model options include "\code{PGM}", "\code{TPGM}", "\code{SPGM}", "\code{GGM}" and "\code{ISM}".}
  \item{graph.type}{graph structure with 3 options:"\code{scale-free}", "\code{hub}", and "\code{lattice}". Default to "\code{scale-free}".}
}
\details{This function will first generate a graph of the specified graph structure; 
 then based on the generated network, it simulates a multivariate data matrix that follows distribution for the Markov Random Fields model specified.}
\value{
A list of two elements:
  \item{B }{\emph{p}x\emph{p} adjacency matrix of the generated graph.}
  \item{X }{\emph{p}x\emph{n} data matrix.}
%% ...
}
\examples{
library(XMRF)

# simulate scale-free network and data of multivariate Poisson for LPGM
sim <- XMRF.Sim(n=100, p=20, model="LPGM", graph.type="scale-free")
hist(sim$X)
plotNet(sim$B)

# simulate hub network and data of multivariate Gaussian for GGM
sim <- XMRF.Sim(n=100, p=20, model="GGM", graph.type="hub")
hist(sim$X)
plotNet(sim$B)

# simulate hub network and data of multivariate bionomial for ISM
sim <- XMRF.Sim(n=100, p=15, model="ISM", graph.type="hub")
hist(sim$X)
plotNet(sim$B)

}
