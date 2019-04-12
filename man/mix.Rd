\name{XTRA}
\alias{mixed}
\alias{mtmixed}
\alias{markov}
\alias{mtgsru}
\title{
Mixed model solver
}
\description{
Function to solve univariate mixed models with or without the usage of omic information. This function allows single-step modeling of replicated observations with marker information available through the usage of a linkage function to connect to a whole-genome regression method. Genomic estimated values can be optionally deregressed (no shrinkage) while fitting the model, which may inflate the variance component associated to this component.
}
\usage{
mixed(y,random=NULL,fixed=NULL,data=NULL,X=list(),alg=emML,maxit=10,Deregress=TRUE,...)
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

where \eqn{y} is the response variable, \eqn{Xb} corresponds to the fixed effect term, \eqn{Zu} correspond to random effect term(s), \eqn{W} is the incidence matrix of terms with omic information and \eqn{a} is omic values by \eqn{a=Mg}, where \eqn{M} is the genotypic matrix and \eqn{g} are marker effects. Here, \eqn{e} is the residual term. An example is provided using the data from the NAM package with: \code{demo(mixedmodel)}

}
\value{
The function wgr returns a list with observed values (\code{obs}), fitted value (\code{hat}), estimated variance components (\code{var}) and regression coefficients (\code{coef}). Regression coefficients would include fixed terms (\code{B}), random terms (\code{U}) and marker effects (\code{G}).
}

\author{
Alencar Xavier
}
