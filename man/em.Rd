\name{WGR2 (EM)}
\alias{emBA}
\alias{emBB}
\alias{emBC}
\alias{emRR}
\alias{emBL}
\alias{emDE}
\alias{emEN}
\title{
Expectation-Maximization WGR
}
\description{
Univariate models to find breeding values through regression fitted via expectation-maximization implemented in C++.
}
\usage{
emRR(y, gen, df = 4, R2 = 0.5)
emBA(y, gen, df = 4, R2 = 0.5)
emBB(y, gen, df = 4, R2 = 0.5, Pi = 0.7)
emBC(y, gen, df = 4, R2 = 0.5, Pi = 0.7)
emBL(y, gen, R2 = 0.5, alpha = 0.02)
emEN(y, gen, R2 = 0.5, alpha = 0.02)
emDE(y, gen, R2 = 0.5)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is not allowed.
}
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers.
}
  \item{df}{
Hyperprior degrees of freedom of variance components.
}
  \item{R2}{
Expected R2, used to calculate the prior shape as proposed by de los Campos et al. (2013).
}
  \item{Pi}{
Value between 0 and 1. Expected probability pi of having null effect (or 1-Pi if Pi>0.5).
}
  \item{alpha}{
Value between 0 and 1. Intensity of L1 variable selection.
}
}
\details{
These regressions are still under development. The model for the whole-genome regression is as follows:

\deqn{y = mu + Xb + e}

where \eqn{y} is the response variable, \eqn{mu} is the intercept, \eqn{X} is the genotypic matrix, \eqn{b} is the effect of an allele substitution (or regression coefficient) and \eqn{e} is the residual term.

}
\value{
The EM functions returns a list with the intercept (\eqn{mu}), the regression coefficient (\eqn{b}), the fitted value (\eqn{hat}), and the estimated intraclass-correlation (\eqn{h2}).
}

\author{
Alencar Xavier
}
\examples{

data(tpod)

# BLUP
BRR = emRR(y,gen)
cor(y,BRR$hat)

# BayesA
BayesA = emBA(y,gen)
cor(y,BayesA$hat)

# BayesB
BayesB = emBB(y,gen)
cor(y,BayesB$hat)

# BayesC
BayesC = emBC(y,gen)
cor(y,BayesC$hat)

# BayesL 1
BayesL1 = emBL(y,gen)
cor(y,BayesL1$hat)

# BayesL 2
BayesL2 = emDE(y,gen)
cor(y,BayesL2$hat)

# Elastic-Net
ElasticNet = emEN(y,gen)
cor(y,ElasticNet$hat)
}