\name{bWGR-package}
\alias{bWGR-package}
\alias{bWGR}
\docType{package}
\title{
Bayesian Whole-Genome Regression
}
\description{
Whole-genome regression methods on Bayesian framework fitted via EM
    or Gibbs sampling, univariate and multivariate, with optional kernel term and
    sampling techniques.
}
\details{
\tabular{ll}{
Package: \tab bWGR\cr
Type: \tab Package\cr
Version: \tab 1.6.6 \cr
Date: \tab 2020-01-07\cr
License: \tab GPL-3 \cr
}
}
\author{
Alencar Xavier, William Muir, Kyle Kocak, Shizhong Xu, Katy Rainey.

Maintainer: Alencar Xavier <alenxav@gmail.com>
}
\examples{
    \dontrun{
data(tpod)
Fit = wgr(y,gen)
cor(y,Fit$hat)
  }
}
