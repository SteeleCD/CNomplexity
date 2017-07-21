\name{genomeDoubling} 
\alias{genomeDoubling} 
\title{CNomplexity::genomeDoubling} 
\usage{genomeDoubling(segFile,armFile=NULL,nReps=1000,
                       scoreFUN=ABSOLUTEscore,doSim=TRUE,
                       withCNb=TRUE,sampleCol = 1,aCol = 9,
                       bCol = 10,startCol = 4,endCol = 5,
                       chromCol = 3,totCol = 11,head=TRUE,
                       diag=FALSE,sepP=FALSE,doParallel=FALSE,
                       nCores = NULL) }
\arguments{
\item{segFile}{Segment file.}
\item{armFile}{Cytoband file. If \code{NULL}, default will be loaded.}
\item{nReps=1000}{Number of simulations.}
\item{scoreFUN=ABSOLUTEscore}{Score to get p value.}
\item{doSim=TRUE}{Whether to simulate.}
\item{withCNb=TRUE}{Does seg file have b CN?.}
\item{sampleCol}{Seg sample column.}
\item{aCol}{Seg a CN column.}
\item{bCol}{Seg b CN column.}
\item{startCol}{Seg start column.}
\item{endCol}{Seg end column.}
\item{chromCol}{Seg chrom column.}
\item{totCol}{Seg total CN column.}
\item{head=TRUE}{Does seg file have headers?.}
\item{diag=FALSE}{Whether to output diagnosis objects..}
\item{sepP=FALSE}{Vars separate for each allele?.}
\item{doParallel=FALSE}{Run in parallel?.}
\item{nCores}{Number of cores.}
}
\value{
Vector of p-values for genome doubling in each sample.
} 
\description{
Identify samples that are genome doubled.
}
\details{
  Implements tests for genome doubling from Dewhurst et al. (2013) and Carter et al. (2012).
}

