\name{XTRA 1}
\alias{mixed}
\alias{mtmixed}
\alias{mtgsru}

\alias{mm}
\alias{NNS}
\alias{GSFLM}
\alias{GSRR}
\alias{GS2EIGEN}
\alias{NNSEARCH}


\title{
Mixed model solver
}
\description{
Function to solve univariate mixed models with or without the usage of omic information. This function allows single-step modeling of replicated observations with marker information available through the usage of a linkage function to connect to a whole-genome regression method. Genomic estimated values can be optionally deregressed (no shrinkage) while fitting the model.
}
\usage{
mixed(y,random=NULL,fixed=NULL,data=NULL,X=list(),
      alg=emML,maxit=10,Deregress=FALSE,...)
}
\arguments{
  \item{y}{
Response variable from the data frame containg the dataset.
}
  \item{random}{
Formula. Right-hand side formula of random effects.
}
  \item{fixed}{
Formula. Right-hand side formula of fixed effects.
}
  \item{data}{
Data frame containing the response variable, random and fixed terms.
}
  \item{X}{
List of omic incidence matrix. Row names of these matrices connect the omic information to the levels of the indicated random terms (eg. \code{X=list("ID"=gen)}).
}
  \item{alg}{
Function. Whole-genome regression algorithm utilized to solve link functions. These include MCMC (\code{wgr}, \code{BayesB}, etc) and EM (\code{emEN}, \code{emDE}, etc) algorithms. By default, it runs maximum likelihood \code{emML}.
}
\item{maxit}{
Integer. Maximum number of iterations.
}
\item{Deregress}{
Logical. Deregress (unshrink) coefficients while fitting the model?
} 
\item{...}{
Additional arguments to be passed to the whole-genome regression algorithms especified on \code{alg}.
}
 
}
\details{
The model for the whole-genome regression is as follows:

\deqn{y = Xb + Zu + Wa + e}

where \eqn{y} is the response variable, \eqn{Xb} corresponds to the fixed effect term, \eqn{Zu} corresponds to one or more random effect terms, \eqn{W} is the incidence matrix of terms with omic information and \eqn{a} is omic values by \eqn{a=Mg}, where \eqn{M} is the genotypic matrix and \eqn{g} are marker effects. Here, \eqn{e} is the residual term. An example is provided using the data from the NAM package with: \code{demo(mixedmodel)}.

Alterinative (and updated) implementations have similar syntax:
    
    01) \code{mm(y,random=NULL,fixed=NULL,data=NULL,
              M=NULL,bin=FALSE,AM=NULL,it=10,verb=TRUE,
              FLM=TRUE,wgtM=TRUE,cntM=TRUE,nPc=3)}
              
    02) \code{mtmixed = function(resp, random=NULL, fixed=NULL,
              data, X=list(), maxit=10, init=10, regVC=FALSE)}

}
\value{
The function wgr returns a list with Fitness values (\code{Fitness}) containing observation \code{obs}, fitted values \code{hat}, residuals \code{res}, and fitted values by model term \code{fits}; Estimated variance components (\code{VarComp}) containing the variance components per se (\code{VarComponents}) and variance explained by each model term (\code{VarExplained}), regression coefficients by model term (\code{Coefficients}), and the effects of structured terms (\code{Structure}) containing the marker effects of each model term where markers were provided.
}

\references{
Xavier, A. (2019). Efficient Estimation of Marker Effects in Plant Breeding. G3: Genes, Genomes, Genetics, DOI: 10.1534/g3.119.400728      
}
\author{
Alencar Xavier
}
 \examples{
 \dontrun{
 demo(mixedmodel)
 }
}
