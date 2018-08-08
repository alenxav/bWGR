\name{bWGR-package}
\alias{bWGR-package}
\alias{bWGR}
\docType{package}
\title{
Bayesian Whole-Genome Regression
}
\description{
Whole-genome regression methods on Bayesian framework fitted via EM
or Gibbs sampling, with optional sampling techniques and kernel term.
}
\details{
\tabular{ll}{
Package: \tab bWGR\cr
Type: \tab Package\cr
Version: \tab 1.5.4 \cr
Date: \tab 2018-08-08\cr
License: \tab GPL-3 \cr
}
}
\author{
Alencar Xavier, Shizhong Xu, William Muir, Katy Rainey.

Maintainer: Alencar Xavier <alenxav@gmail.com>
}
\examples{
# Fit BayesB
data(tpod)
Fit = BayesB(y,gen)
cor(y,Fit$hat)
}
