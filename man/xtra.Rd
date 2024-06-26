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
\alias{AccByC}
\alias{EigenArcZ}
\alias{EigenGauZ}
\alias{K2X}
\alias{SimY}
\alias{SimZ}
\alias{SimGC}
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

\code{GRM(X,Code012=FALSE)} # Creates additive kinship matrix VanRaden 2008

\code{SPC(y,blk,row,col,rN=3,cN=1)} # Spatial covariate

\code{SPM(blk,row,col,rN=3,cN=1)} # Spatial design matrix

\code{SibZ(id,p1,p2)} # Pedigree design matrix compatible to regression methods

\code{Hmat(ped,gen=NULL)} # Kinship combining pedigree and genomics
  
\code{EigenGRM(X, centralizeZ = TRUE, cores = 1)} # GRM using Eigen library

\code{EigenARC(X, centralizeZ = TRUE, cores = 1)} # ArcCosine kernel

\code{EigenGAU(X, phi = 1.0, cores = 1)} # Gaussian kernel using Eigen library

\code{EigenCNT(X, cores = 1)} # Center SNPs without missing Eigen library

\code{EigenEVD(A, cores = 1)} # Eigendecomposition from Eigen library

\code{EigenBDCSVD(X, cores = 1)} # BDC single value composition from Eigen

\code{EigenJacobiSVD(X, cores = 1)} # Jacobi single value composition from Eigen

\code{EigenAcc(X1, X2, h2 = 0.5, cores = 1)} # Deterministic accuracy X1 -> X2 via V

\code{AccByC(X1, X2, h2 = 0.5, cores = 1)} # Deterministic accuracy X1 -> X2 via C

\code{EigenArcZ(Zfndr, Zsamp, cores = 1)} # Reduced rank ArcCos kernel PCs with founder rotation

\code{EigenGauZ(Zfndr, Zsamp, phi=1, cores = 1)} # Reduced rank Gaussian kernel PCs with founder rotation 

\code{K2X(K, MinEV = 1e-8, cores = 1)} # Reparametrize kernel to PCs to run regression models 

\code{SimY(Z,k=5,h2=0.5,GC=0.5,seed=123,unbalanced=FALSE,PercMiss=0,BlkMiss=FALSE)} # Simulate phenotypes 

\code{SimZ(ind=500,snp=500,chr=2,F2=TRUE,rec=0.01)} # Simulate genome

\code{SimGC(k=50,...)} # Simulate genetic correlation matrix

\code{MvSimY(Ufndr,Zfndr,Zsamp,GxY,GxL,H2plot,nLoc=20,Seed=123)} # Simulate phenotypes given founders

}
\author{
Alencar Xavier
}
