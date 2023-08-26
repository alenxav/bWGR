\name{XTRA 2}
\alias{emGWA}
\alias{SPC}
\alias{SPM}
\alias{markov}
\alias{IMP}
\alias{CNT}
\alias{GAU}
\alias{GRM}
\alias{Hmat}
\alias{SibZ}
\alias{EigenGAU}
\alias{EigenGRM}
\alias{EigenCNT}
\alias{EigenARC}
\alias{EigenEVD}
\alias{EigenBDCSVD}
\alias{EigenJacobiSVD}
\alias{EigenAcc}
\alias{EigenArcZ}
\alias{EigenGauZ}
\alias{K2X}
\alias{SimY}
\alias{MvSimY}

\title{
Additional tools
}
\description{
Complementary functions that may help with handling parameters and routine operations.
}

\details{

\code{emGWA(y,gen)} # Simple MLM for association analysis

\code{markov(gen,chr=NULL)} # Markovian imputation of genotypes coded as 012

\code{IMP(X)} # Imputes genotypes with SNP expectation (column average)

\code{CNT(X)} # Recodes SNPs by centralizing columns in a matrix

\code{GAU(X)} # Creates Gaussian kernel as exp(-Dist2/mean(Dist2))

\code{GRM(X,Code012=FALSE)} # Creates additive kinship matrix VanRaden (2008)

\code{SPC(y,blk,row,col,rN=3,cN=1)} # Spatial covariate

\code{SPM(blk,row,col,rN=3,cN=1)} # Spatial design matrix

\code{SibZ(id,p1,p2)} # Pedigree design matrix compatible to regression methods

\code{Hmat(ped,gen=NULL)} # Kinship combining pedigree and genomics
  
\code{EigenGRM(X, centralizeZ = TRUE, cores = 2)} # GRM using Eigen library

\code{EigenARC(X, centralizeZ = TRUE, cores = 2)} # ArcCosine kernel

\code{EigenGAU(X, phi = 1.0, cores = 2)} # Gaussian kernel using Eigen library

\code{EigenCNT(X, cores = 2)} # Center SNPs without missing Eigen library

\code{EigenEVD(A, cores = 2)} # Eigendecomposition from Eigen library

\code{EigenBDCSVD(X, cores = 2)} # BDC single value composition from Eigen

\code{EigenJacobiSVD(X, cores = 2)} # Jacobi single value composition from Eigen

\code{EigenAcc(X1, X2, h2 = 0.5, cores = 2)} # Deterministic accuracy X1 -> X2

\code{EigenArcZ(Zfndr, Zsamp, cores = 2)} # Reduced rank ArcCos kernel PCs with founder rotation

\code{EigenGauZ(Zfndr, Zsamp, phi=1, cores = 2)} # Reduced rank Gaussian kernel PCs with founder rotation 

\code{K2X(K)} # Reparametrize kernel to PCs to run regression models 

\code{SimY(Z,k=5,h2=0.5,GC=0.5,seed=123,unbalanced=FALSE)} # Simulate Y for multiple environments

\code{MvSimY(Ufndr,Zfndr,Zsamp,GxY,GxL,H2plot,nLoc=20,Seed=123)} # Simulate traits given founders



}
\author{
Alencar Xavier
}
