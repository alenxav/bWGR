\name{WGR2 (EM)}
\alias{emBA}
\alias{emBB}
\alias{emBC}
\alias{emRR}
\alias{emBL}
\alias{emDE}
\alias{emEN}
\alias{emCV}
\alias{SPC}
\alias{IMP}
\alias{CNT}
\alias{GAU}
\title{
Expectation-Maximization WGR
}
\description{
Univariate models to find breeding values through regression fitted via expectation-maximization implemented in C++.
}
\usage{
emRR(y, gen, df = 10, R2 = 0.5)
emBA(y, gen, df = 10, R2 = 0.5)
emBB(y, gen, df = 10, R2 = 0.5, Pi = 0.75)
emBC(y, gen, df = 10, R2 = 0.5, Pi = 0.75)
emBL(y, gen, R2 = 0.5, alpha = 0.02)
emEN(y, gen, R2 = 0.5, alpha = 0.02)
emDE(y, gen, R2 = 0.5)
emCV(y, gen, k = 5, n = 5)
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
  \item{k}{
Integer. Folding of a k-fold cross-validation.
}
  \item{n}{
Integer. Number of cross-validation to perform.
}
}
\details{
The model for the whole-genome regression is as follows:

\deqn{y = mu + Xb + e}

where \eqn{y} is the response variable, \eqn{mu} is the intercept, \eqn{X} is the genotypic matrix, \eqn{b} is the effect of an allele substitution (or regression coefficient) and \eqn{e} is the residual term. A k-fold cross-validation for model evaluation is provided by \eqn{emCV}.

}
\value{
The EM functions returns a list with the intercept (\eqn{mu}), the regression coefficient (\eqn{b}), the fitted value (\eqn{hat}), and the estimated intraclass-correlation (\eqn{h2}). 

The function emCV returns the predictive ability of each model, that is, the correlation between the predicted and observed values from \eqn{k}-fold cross-validations repeated \eqn{n} times.
}

\author{
Alencar Xavier
}
\examples{

data(tpod)
emCV(y,gen,3,3)

}