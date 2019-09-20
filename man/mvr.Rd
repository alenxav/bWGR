\name{WGR3 (MV)}
\alias{mkr}
\alias{mrr}
\alias{mrrR}
\alias{mkr2X}
\alias{mrr2X}
\alias{mrrFast}
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
    
    \deqn{Y = Mu + UB + E}
  
  where \eqn{Y} is a matrix of response variables, \eqn{Mu} represents the intercepts, \eqn{U} is the matrix of Eigenvector of K, \eqn{b} is a vector of regression coefficients and \eqn{E} is the residual matrix.
  
  The model for the ridge regression (mrr) is as follows:
    
    \deqn{Y = Mu + XB + E}
  
  where \eqn{Y} is a matrix of response variables, \eqn{Mu} represents the intercepts, \eqn{X} is the matrix of genotypic information, \eqn{B} is the matrix of marker effects,  and \eqn{E} is the residual matrix.
  
  Algorithm: Residuals are assumed to be independent among traits. Regression coefficients are solved via a multivaraite adaptation of Gauss-Seidel Residual Update. Variance and covariance components are solved with an efficient variation of EM-REML.
  
  Other related implementations:
    
    01) \code{mkr2X(Y,K1,K2):}{ Solves multi-trait kernel regressions with two random effects.}

    02) \code{mrr2X(Y,X1,X2):}{ Solves multi-trait ridge regressions with two random effects.}

    03) \code{mrrR(Y,X):}{ Variation of \code{mrr} that assumes correlated residuals. }

}
\value{
  
  Returns a list with the random effect covariances (\code{Vb}), residual variances (\code{Ve}), genetic correlations (\code{GC}), matrix with marker effects (\code{b}) or eigenvector effects (if \code{mkr}), intercepts (\code{mu}), heritabilities (\code{h2}), and a matrix with fitted values (\code{hat}).
  
}

\author{
  Alencar Xavier
}
\examples{
  \dontrun{
    
    # Load data and compute G matrix
    data(tpod)
    gen = CNT(gen)
    K = tcrossprod(gen)
    K = K/mean(diag(K))
    
    # Phenotypes: 3 traits correlated r=0.5
    G0 = 0.5+diag(0.5,3) 
    G = kronecker(G0,K)
    diag(G)=diag(G)+0.001
    L = chol(G)
    TBV = crossprod(L,rnorm(196*3))
    Y = rnorm(196*3,10+TBV,sd(TBV))
    Phe = matrix(Y,ncol=3)
    TBV = matrix(TBV,ncol=3)
    
    # Fit kernel and regression models
    test1 = mkr(Phe,K)
    test2 = mrr(Phe,gen)
    
    # Genetic correlation
    test1$GC
    test2$GC
    
    # Heritabilies
    test1$h2
    test2$h2
    
    # Goodness of fit
    diag(cor(TBV,test1$hat))
    diag(cor(TBV,test2$hat))
    
    
  }
}
