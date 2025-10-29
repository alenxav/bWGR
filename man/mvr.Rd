\name{WGR3 (MV)}
\alias{mkr}
\alias{mrr}
\alias{mrr_float}
\alias{mkr2X}
\alias{mrr2X}
\alias{MRR3}
\alias{MRR3F}
\alias{MLM}
\alias{SEM}
\alias{mrr_svd}
\alias{MEGA}
\alias{GSEM}
\alias{UVBETA}
\alias{solver1x}
\alias{solver2x}
\alias{FUVBETA}
\alias{solver1xF}
\alias{solver2xF}
\alias{FUVBETA}
\alias{XFUVBETA}
\alias{ZFUVBETA}
\alias{XSEMF}
\alias{YSEMF}
\alias{ZSEMF}
\alias{PEGS}
\alias{mwgr}
\alias{ClusterBlup}
\alias{IncMatrix}
\alias{EM_recluster}
\alias{Get_Cluster_Corr}

\title{
  Multivariate Regression
}
\description{
  Multivariate model to find breeding values.
}
\usage{
  mkr(Y,K,...)
  mrr(Y,X,...)
  mrr_float(Y,X,...)
  mwgr(Y,X,solver=MRR3F,...)
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
  \item{...}{
    Arguments to pass to solver (e.g., MRR3). See \code{args(MRR3)}.
  } 
  \item{solver}{
    Which multivariate solver to use for mwgr wrapper.
  } 
  
}
\details{

 Multivariate solver algorithms are described in Xavier and Habier (2022) and Xavier et al. (2024).
 
 The mwgr function is a wrapper that adds row and column names to the output of multivariate solvers.
 
 The model for the ridge regression (mrr) is as follows:
    
    \deqn{Y = Mu + XB + E}
  
  where \eqn{Y} is a matrix of response variables, \eqn{Mu} represents the intercepts, \eqn{X} is the matrix of genotypic information, \eqn{B} is the matrix of marker effects,  and \eqn{E} is the residual matrix.
  
  The model for the kernel regression (mkr) is as follows:
    
    \deqn{Y = Mu + UB + E}
  
  where \eqn{Y} is a matrix of response variables, \eqn{Mu} represents the intercepts, \eqn{U} is the matrix of Eigenvector of K, \eqn{b} is a vector of regression coefficients and \eqn{E} is the residual matrix.
  
   Algorithm: Residuals are assumed to be independent among traits. Regression coefficients are solved via a multivaraite adaptation of Gauss-Seidel Residual Update. Since version 2.0, the solver of \code{mrr} is based on the Randomized Gauss-Seidel algorithm. Variance and covariance components are solved with an EM-REML like approach proposed by Schaeffer called Pseudo-Expectation.
  
  Other related implementations:
    
    01) \code{mkr2X(Y,K1,K2):} Solves multi-trait kernel regressions with two random effects.

    02) \code{mrr2X(Y,X1,X2):} Solves multi-trait ridge regressions with two random effects.

    03) \code{MRR3(Y,X,...):} Extension of mrr with additional parameters.

    04) \code{MRR3F(Y,X,...):} MRR3 running on float.

    05) \code{mrr_svd(Y,W):} Solves mrr through the principal components of parameters.

    06) \code{MLM(Y,X,Z,maxit=500,logtol=-8,cores=1):} Multivariate model with fixed effects.

    07) \code{SEM(Y,Z,...):} Fits a MegaSEM with both shared- and trait-specific terms.

    08) \code{MEGA(Y,X,npc=-1):} Toy implementation of MegaLMM, imputing missing with GEBVs.

    09) \code{GSEM(Y,X,npc=-1):} Toy C++ implementaiton of MegaSEM, jointly fits FA and XB.

    10) \code{ZSEMF(Y,X,npc=0):} Full-rank MegaSEM, float precision.

    11) \code{YSEMF(Y,X,npc=-1):} Reduced-rank MegaSEM, float, two-steps approach.

    12) \code{XSEMF(Y,X,npc=0):} Full-rank MegaSEM, h2 fixed at 0.5, float precision.
    
    13) \code{PEGS(Y,X,...):} Light PEGS implementation, float precision.
    
In GSEM, XSEMF and MEGA, 'npc' means number of latent spaces if input is above zero, otherwise, 0 means all and -1 means \code{2*sqrt(ncol(Y))}.

Additional functions (ClusterBlup, IncMatrix, EM_recluster, Get_Cluster_Corr) may help with downstream analysis.
    
}
\value{
  
  Returns a list with the random effect covariances (\code{Vb}), residual variances (\code{Ve}), genetic correlations (\code{GC}), matrix with marker effects (\code{b}) or eigenvector effects (if \code{mkr}), intercepts (\code{mu}), heritabilities (\code{h2}), and a matrix with fitted values (\code{hat}).
  
   NOTE: Numeric stability is a serious concern with multivariate models with large number of response variables, as the covariance matrix is often not invesible. If output is filled with NAs, try using \code{MRR3} and play with some parameters. For example, one may try adding priors to stabilize variances, e.g., \code{fit=MRR3(Y,X,df0=20)}. 
  
  
  
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
    
    # try: demo(multivariates)

}

\references{
Xavier, A and Habier, D. (2022). A new approach fits multivariate genomic prediction models efficiently. GSE, DOI: 10.1186/s12711-022-00730-w 

Xavier, A et al. (2024). Megavariate Methods Capture Complex Genotype-by-Environment Interactions. Genetics, DOI: 10.1093/genetics/iyae179
}
