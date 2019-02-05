\name{WGR3 (MVR)}
\alias{mkr}
\alias{mrr}
\alias{Hmat}
\title{
Multi-variate Regression
}
\description{
Multivariate model to find breeding values.
}
\usage{
 mkr(Y,K=NULL,eK=NULL,it=500,bu=200,th=3,df=5,R2=0.5,EigT=0.05,verb=FALSE)
 mrr(Y,X,Choleski=FALSE)    
}
\arguments{
  \item{Y}{
Numeric matrix of observations (\eqn{n,k}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{X}{
Numeric matrix containing the genotyping matrix.
} 
  \item{Choleski}{
Use Choleski decomposition to invert genetic covariance.
} 
  \item{K}{
Numeric matrix containing the genotypic relationship matrix. A matrix with \eqn{n} rows and columns.
}
  \item{eK}{
Output of eigen. Eigendecomposition of K. If eK is specified, the is no need for the argument K.
}
  \item{it}{
Integer. Number of iterations or samples to be generated.
}
  \item{bu}{
Integer. Burn-in, the number of iterations or samples to be discarted.
}
  \item{th}{
Integer. Thinning parameter, used to save memory by storing only one every 'th' samples.
}
  \item{df}{
Prior degrees of freedom for covariance components.
}
  \item{R2}{
Expected R2, used to calculate the prior shape as proposed by de los Campos et al. (2013).
}
  \item{EigT}{
Null or numeric. If provided, the model uses just Eigenpairs with Eigenvalues above the specified theshold.
}
  \item{verb}{
Logical. If verbose is TRUE, function displays MCMC progress bar.
}

}
\details{
The model for the kernel regression is as follows:

\eqn{Y = mu + Z(UB) + E},

where \eqn{Y} is a matrix of response variables, \eqn{mu} represents the intercepts, \eqn{Z} is the design matrix, \eqn{U} is the matrix of Eigenvector of K, \eqn{b} is a vector of regression coefficients and \eqn{E} is the residual matrix. Variance components are sampled from a inverse Wishart distribution (Sorensen and Gianola 2002). Regression coefficients are solved with an adaptation of the algorithm proposed by de los Campos et al. (2010).

}
\value{
The function mkr returns a list with the random effect covariance matrix (\eqn{VA}), residual covariance matrix (\eqn{VE}) and a matrix with breeding values (\eqn{BV}).
}

\references{

de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345.

de los Campos, G., Gianola, D., Rosa, G. J., Weigel, K. A., & Crossa, J. (2010). Semi-parametric genomic-enabled prediction of genetic values using reproducing kernel Hilbert spaces methods. Genetics Research, 92(04), 295-308.

Sorensen D., and Gianola D. (2002) Likelihood, Bayesian, and MCMC methods in quantitative genetics. Springer.

}
\author{
Alencar Xavier
}
\examples{
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
test2 = mrr(Phe,gen,TRUE)

# Genetic correlation
cov2cor(test1$VA)

# Environmental correlation
cov2cor(test1$VE)
 
# Heritabilies
diag(test1$VA/(test1$VA+test1$VE))
test2$h2
 
# Goodness of fit
diag(cor(Phe,test1$BV))
diag(cor(Phe,test2$hat))
 
}
