\name{bWGR-package}
\alias{bWGR-package}
\alias{bWGR}
\docType{package}
\title{
Bayesian Whole-Genome Regression
}
\description{
Whole-genome regression methods on Bayesian framework fitted via EM
    or Gibbs sampling, with optional kernel term and sampling techniques as described
    by Xavier (2017) <doi:10.1186/s12859-017-1582-3>.
}
\details{
\tabular{ll}{
Package: \tab bWGR\cr
Type: \tab Package\cr
Version: \tab 1.6.5 \cr
Date: \tab 2019-05-20\cr
License: \tab GPL-3 \cr
}
}
\author{
Alencar Xavier, William Muir, Kyle Kocak, Shizhong Xu, Katy Rainey.

Maintainer: Alencar Xavier <alenxav@gmail.com>
}
\examples{
    \dontrun{
# Fit Bayesian Ridge Regression
data(tpod)
Fit = wgr(y,gen)
cor(y,Fit$hat)
  }
}
