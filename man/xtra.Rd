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
}
\author{
Alencar Xavier
}
