\name{WGR1 (MCMC)}
\alias{wgr}
\alias{KMUP}
\alias{KMUP2}
\alias{GSEN}
\alias{BayesA}
\alias{BayesA2}
\alias{BayesB}
\alias{BayesB2}
\alias{BayesC}
\alias{BayesCpi}
\alias{BayesDpi}
\alias{BayesL}
\alias{BayesRR}
\alias{BayesRR2}
\alias{mcmcCV}
\title{
Whole-genome Regression
}
\description{
Univariate model to find breeding values through regression with optional resampling techniques (SBMC) and polygenic term (Kernel). Additional standalone functions written entirely in C++ are also implemented.
}
\usage{
wgr(y,X,it=1500,bi=500,th=1,bag=1,rp=FALSE,iv=FALSE,de=FALSE,
    pi=0,df=5,R2=0.5,eigK=NULL,VarK=0.95,verb=FALSE)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{X}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers.
}
  \item{it}{
Integer. Number of iterations or samples to be generated.
}
  \item{bi}{
Integer. Burn-in, the number of iterations or samples to be discarted.
}
  \item{th}{
Integer. Thinning parameter, used to save memory by storing only one every 'th' samples.
}
  \item{bag}{
If different than one, it indicates the proportion of data to be subsampled in each Markov chain. For datasets with moderate number of observations, values of bag from 0.30 to 0.60 may speed up computation without losses in predicion properties. This argument enable users to enhance MCMC through Subsampling bootstrap Markov chain (SBMC).
}
  \item{rp}{
Logical. Use replacement for bootstrap samples when bag is different than one.
}
  \item{iv}{
Logical. Assign markers independent variance, a T prior from a mixture of normals. If true, turns the default model BLUP into BayesA.
}
  \item{de}{
Logical. Assign markers independent variance through double-exponential prior. If true, turns the default model BLUP into Bayesian LASSO. This argument overides iv.
}
  \item{pi}{
Value between 0 and 1. If greater than zero it activates variable selection, where markers have expected probability pi of having null effect.
}
  \item{df}{
Prior degrees of freedom of variance components.
}
  \item{R2}{
Expected R2, used to calculate the prior shape.
}
  \item{eigK}{
Output of function 'eigen'. Spectral decomposition of the kernel used to compute the polygenic term.
}
  \item{VarK}{
Numeric between 0 and 1. For reduction of dimensionality. Indicates the proportion of variance explained by Eigenpairs used to fit the polygenic term.
}
  \item{verb}{
Logical. If verbose is TRUE, function displays MCMC progress bar.
}
}
\details{
The model for the whole-genome regression is as follows:

\deqn{y = mu + Xb + u + e}

where \eqn{y} is the response variable, \eqn{mu} is the intercept, \eqn{X} is the genotypic matrix, \eqn{b} is the regression coefficient or effect of an allele substitution, with \eqn{d} probability of being included into the model, \eqn{u} is the polygenic term if a kernel is used, and \eqn{e} is the residual term.

Users can obtain four WGR methods out of this function: BRR (pi=0,iv=F), BayesA (pi=0,iv=T), BayesB (pi=0.95,iv=T), BayesC (pi=0.95,iv=F) and Bayesian LASSO or BayesL (pi=0,de=T). Theoretical basis of each model is described by de los Campos et al. (2013).

Gibbs sampler that updates regression coefficients is adapted from GSRU algorithm (Legarra and Misztal 2008). The variable selection of functions \eqn{wgr}, \eqn{BayesB} and \eqn{BayesC} works through the unconditional prior algorithm proposed by Kuo and Mallick (1998), whereas \eqn{BayesCpi} and \eqn{BayesDpi} are based on Metropolis-Hastings. Prior shape estimates are computed as Sb = R2*df*var(y)/MSx and Se = (1-R2)*df*var(y), with an exception for \eqn{BayesC} and \eqn{BayesCpi} where the prior shape is Sb = R2*df*var(y)/MSx/(1-pi). The polygenic term is solved by Bayesian algorithm of reproducing kernel Hilbert Spaces proposed by de los Campos et al. (2010).


In addition to \code{wgr}, standalone C++ functions available include:

01) \code{BayesA(y,X,it=1500,bi=500,df=5,R2=0.5)}

02) \code{BayesB(y,X,it=1500,bi=500,pi=0.95,df=5,R2=0.5)}

03) \code{BayesC(y,X,it=1500,bi=500,pi=0.95,df=5,R2=0.5)}

04) \code{BayesCpi(y,X,it=1500,bi=500,df=5,R2=0.5)}
    
05) \code{BayesDpi(y,X,it=1500,bi=500,df=5,R2=0.5)}
    
06) \code{BayesL(y,X,it=1500,bi=500,df=5,R2=0.5)}

07) \code{BayesRR(y,X,it=1500,bi=500,df=5,R2=0.5)}

The implementations that support two random effects include:

08) \code{BayesA2(y,X1,X2,it=1500,bi=500,df=5,R2=0.5)}

09) \code{BayesB2(y,X1,X2,it=1500,bi=500,pi=0.95,df=5,R2=0.5)}

10) \code{BayesRR2(y,X1,X2,it=1500,bi=500,df=5,R2=0.5)}

And the cross-validation for the C++ implementations, with arguments analogous to the emCV function.
    
Cross-validation: \code{mcmcCV(y,gen,k=5,n=5,it=1500,bi=500,pi=0.95,df=5,R2=0.5,avg=TRUE,llo=NULL,tbv=NULL)}    
    
}
\value{
The function wgr returns a list with expected value from the marker effect (\eqn{b}), probability of marker being in the model (\eqn{d}), regression coefficient (\eqn{g}), variance of each marker (\eqn{Vb}), the intercept (\eqn{mu}), the polygene (\eqn{u}) and polygenic variance (\eqn{Vk}), residual variance (\eqn{Ve}) and the fitted value (\eqn{hat}).
}
\references{

de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345.

de los Campos, G., Gianola, D., Rosa, G. J., Weigel, K. A., & Crossa, J. (2010). Semi-parametric genomic-enabled prediction of genetic values using reproducing kernel Hilbert spaces methods. Genetics Research, 92(04), 295-308.

Kuo, L., & Mallick, B. (1998). Variable selection for regression models. Sankhya: The Indian Journal of Statistics, Series B, 65-81.

Legarra, A., & Misztal, I. (2008). Technical note: Computing strategies in genome-wide selection. Journal of dairy science, 91(1), 360-366.

}
\author{
Alencar Xavier
}
\examples{
    \dontrun{
        
# Load data
data(tpod)

# BLUP
fit_BRR = wgr(y,gen,iv=FALSE,pi=0)
cor(y,fit_BRR$hat)

# BayesA
fit_BayesA = wgr(y,gen,iv=TRUE,pi=0)
cor(y,fit_BayesA$hat)

# BayesB
fit_BayesB = wgr(y,gen,iv=TRUE,pi=.95)
cor(y,fit_BayesB$hat)

# BayesC
fit_BayesC = wgr(y,gen,iv=FALSE,pi=.95)
cor(y,fit_BayesC$hat)

# BayesCpi
fit_BayesCpi = BayesCpi(y,gen)
cor(y,fit_BayesCpi$hat)        

# BayesCpi
fit_BayesDpi = BayesDpi(y,gen)
cor(y,fit_BayesDpi$hat)        
        
# BayesL
fit_BayesL = wgr(y,gen,de=TRUE)
cor(y,fit_BayesL$hat)

# Bagging BLUP
fit_Bag = wgr(y,gen,bag=0.5)
cor(y,fit_Bag$hat)
    
    }
}
