\name{WGR2 (EM)}
\alias{emBA}
\alias{emBB}
\alias{emBC}
\alias{emRR}
\alias{emBL}
\alias{emDE}
\alias{emEN}
\alias{emML}
\alias{emML2}
\alias{emCV}
\alias{emBCpi}
\alias{lasso}

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
emBCpi(y, gen, df = 10, R2 = 0.5, Pi = 0.75)
emBL(y, gen, R2 = 0.5, alpha = 0.02)
emEN(y, gen, R2 = 0.5, alpha = 0.02)
emDE(y, gen, R2 = 0.5)
emML(y, gen, D = NULL)
lasso(y, gen)

emCV(y, gen, k = 5, n = 5, Pi = 0.75, alpha = 0.02,
     df = 10, R2 = 0.5, avg=TRUE, llo=NULL, tbv=NULL, ReturnGebv = FALSE)
}
\arguments{
  \item{y}{
Numeric vector of response variable (\eqn{n}). \code{NA} is not allowed.
}
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and \eqn{m} columns of molecular markers.
}
  \item{df}{
Hyperprior degrees of freedom of variance components.
}
  \item{R2}{
Expected R2, used to calculate the prior shape (de los Campos et al. 2013).
}
  \item{Pi}{
Value between 0 and 1. Expected probability pi of having null effect (or 1-Pi if Pi>0.5).
}
  \item{alpha}{
Value between 0 and 1. Intensity of L1 variable selection.
}
  \item{D}{
NULL or numeric vector with length p. Vector of weights for markers.
}
  \item{k}{
Integer. Folding of a k-fold cross-validation.
}
  \item{n}{
Integer. Number of cross-validation to perform.
}
  \item{avg}{
Logical. Return average across CV, or correlations within CV.
}
  \item{llo}{
NULL or a vector (numeric or factor) with the same length as y. If provided, the cross-validations are performed as Leave a Level Out (LLO). This argument allows the user to predefine the splits. This argument overrides \code{k} and \code{n}.
}
  \item{tbv}{
NULL or numeric vector of 'true breeding values' (\eqn{n}) to use to compare cross-validations to. If NULL, the cross-validations will have the phenotypes as prediction target.
}
  \item{ReturnGebv}{
Logical. If TRUE, it returns a list with the average marker values and fitted values across all cross-validations, in addition to the regular output.
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
     \dontrun{

data(tpod)
emCV(y,gen,3,3)
          
 }
}
