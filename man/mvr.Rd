\name{WGR3 (MV)}
\alias{mkr}
\alias{mkr_old}
\alias{mkr2X}
\alias{mrr}
\alias{mrrV2}
\alias{mrrV3}
\alias{mrr2X}
\alias{Hmat}
\alias{SibZ}
\title{
Multivariate Regression
}
\description{
Multivariate model to find breeding values.
}
\usage{
 mkr(Y,K)
 mrr(Y,X)
}
\arguments{
  \item{Y}{
Numeric matrix of observations x trait. \code{NA} is allowed.
}
  \item{K}{
Numeric matrix containing the relationship matrix.
}
 \item{X}{
Numeric matrix containing the genotyping matrix.
} 

}
\details{

The model for the kernel regression (mkr) is as follows:

\eqn{Y = Mu + UB + E},

where \eqn{Y} is a matrix of response variables, \eqn{Mu} represents the intercepts, \eqn{U} is the matrix of Eigenvector of K, \eqn{b} is a vector of regression coefficients and \eqn{E} is the residual matrix.

The model for the ridge regression (mrr) is as follows:

\eqn{Y = Mu + XB + E},

where \eqn{Y} is a matrix of response variables, \eqn{Mu} represents the intercepts, \eqn{X} is the matrix of genotypic information, \eqn{B} is the matrix of marker effects,  and \eqn{E} is the residual matrix. If \code{mrr} returns an error, an alternative implementation (\code{mrrV2}) in available. 

 Algorithm: Residuals are assumed to be independent among traits. Regression coefficients are solved via a multivaraite adaptation of Gauss-Seidel Residual Update. Covariance components are solved with an efficient EM-REML.

 Other related implementations:
 
 01) \code{mkr2X(Y,K1,K2):}{ Solves multi-trait kernel regressions with two random effects.}

 02) \code{mrr2X(Y,X1,X2):}{ Solves multi-trait ridge regressions with two random effects. }

 03) \code{mrrV2(Y,X):}{ More stable variation of \code{mrr} which ridges the genetic covariance (add small bias to the diagonal) to assure inversion. }
 
 04) \code{mrrV3(Y,X):}{ Variation of \code{mrr} that assumes correlated residuals. }
 
 05) \code{mkr_old(Y,K=NULL,it=500,bu=200,th=3,df=5,R2=0.5):}{ Old MCMC implementation of \code{mkr}. }
 
 
 
}
\value{

Returns a list with the random effect covariances (\code{Vb}), residual variances (\code{Ve}), genetic correlations (\code{GC}), matrix with marker effects (\code{b}) or eigenvector effects (if \code{mkr}), intercepts (\code{mu}), heritabilities (\code{h2}), and a matrix with fitted values (\code{hat}).
 
}

\author{
Alencar Xavier
}
\examples{
 \dontrun{
 
# G matrix
data(tpod)
G = tcrossprod(gen)
G = G/mean(diag(G))

# Phenotypes
Y1 = rnorm(196,y,.1)
Y2 = rnorm(196,y,.2)
Y3 = rnorm(196,y,.3)
Phe = cbind(Y1,Y2,Y3)

# Fit model
test1 = mkr(Phe,G)
test2 = mrr(Phe,gen)

# Genetic correlation
test1$GC
test2$GC
 
# Heritabilies
test1$h2
test2$h2
 
# Goodness of fit
diag(cor(Phe,test1$hat))
diag(cor(Phe,test2$hat))
  
 }
}
