\name{WGR3 (MV)}
\alias{mkr}
\alias{mrr}
\alias{mkr2X}
\alias{mrr2X}
\alias{MRR3}
\alias{MLM}
\alias{SEM}
\alias{mrr_svd}

\title{
  Multivariate Regression
}
\description{
  Multivariate model to find breeding values.
}
\usage{
  mrr(Y,X)
  mkr(Y,K)
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
  
 Algorithm is described in Xavier and Habier (2022). The model for the ridge regression (mrr) is as follows:
    
    \deqn{Y = Mu + XB + E}
  
  where \eqn{Y} is a matrix of response variables, \eqn{Mu} represents the intercepts, \eqn{X} is the matrix of genotypic information, \eqn{B} is the matrix of marker effects,  and \eqn{E} is the residual matrix.
  
  The model for the kernel regression (mkr) is as follows:
    
    \deqn{Y = Mu + UB + E}
  
  where \eqn{Y} is a matrix of response variables, \eqn{Mu} represents the intercepts, \eqn{U} is the matrix of Eigenvector of K, \eqn{b} is a vector of regression coefficients and \eqn{E} is the residual matrix.
  
   Algorithm: Residuals are assumed to be independent among traits. Regression coefficients are solved via a multivaraite adaptation of Gauss-Seidel Residual Update. Since version 2.0, the solver of \code{mrr} is based on the Randomized Gauss-Seidel algorithm. Variance and covariance components are solved with an EM-REML like approach proposed by Schaeffer called Pseudo-Expectation. S
  
  Other related implementations:
    
    01) \code{mkr2X(Y,K1,K2):}{ Solves multi-trait kernel regressions with two random effects.}

    02) \code{mrr2X(Y,X1,X2):}{ Solves multi-trait ridge regressions with two random effects.}

    03) \code{MRR3(Y,X,...):}{ Similar to mrr, with a couple additional parameters.}

    04) \code{mrr_svd(Y,W):}{ Solves mrr through the principal components of parameters.}

    05) \code{MLM(Y,X,Z,maxit=500,logtol=-8,cores=1):}{ Multivariate model with fixed effects.}

    06) \code{MLMX(Y,X,Zlist):}{ Multivariate with fixed effects and multiple random terms.}

    07) \code{SEM(Y,Z,PCs=3,TOI=NULL,Beta0=NULL):}{ Fits a XFA structural equation model.}

}
\value{
  
  Returns a list with the random effect covariances (\code{Vb}), residual variances (\code{Ve}), genetic correlations (\code{GC}), matrix with marker effects (\code{b}) or eigenvector effects (if \code{mkr}), intercepts (\code{mu}), heritabilities (\code{h2}), and a matrix with fitted values (\code{hat}).
  
}

\author{
  Alencar Xavier, David Habier
}
\examples{
    
    # Load genomic data
    data(tpod)
    X = CNT(gen)
    
    # Simulate phenotyp
    sim = SimY(X)
    Y = sim$Y
    TBV = sim$tbv
    
    # Fit regression model
    test = mrr(Y,X)
    
    # Genetic correlation
    test$GC
    
    # Heritabilies
    test$h2
    
    # Accuracy
    diag(cor(TBV,test$hat))

}

\references{
Xavier, A and Habier, D. (2022). A new approach fits multivariate genomic prediction models efficiently. GSE, DOI: 10.1186/s12711-022-00730-w    
}
