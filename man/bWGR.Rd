\name{bWGR-package}
\alias{bWGR-package}
\alias{bWGR}
\docType{package}
\title{
Bayesian Whole-Genome Regression
}
\description{
Whole-genome regression methods on Bayesian framework fitted via EM
    or Gibbs sampling, single step (<doi:10.1534/g3.119.400728>), univariate and multivariate (<doi:10.1186/s12711-022-00730-w>),
    with optional kernel term and sampling techniques (<doi:10.1186/s12859-017-1582-3>).
}
\details{
\tabular{ll}{
Package: \tab bWGR\cr
Type: \tab Package\cr
Version: \tab 2.2.10 \cr
Date: \tab 2024-09-17\cr
License: \tab GPL-3 \cr
}
}
\author{
Alencar Xavier, William Muir, David Habier, Kyle Kocak, Shizhong Xu, Katy Rainey.
Maintainer: Alencar Xavier <alenxav@gmail.com>
}
\examples{
data(tpod)
Fit = wgr(y,gen)
cor(y,Fit$hat)
}
