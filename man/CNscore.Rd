\name{CNscore} 
\alias{CNscore} 
\title{CNomplexity::CNscore} 
\usage{CNscore(dir,threshLow=0.275,
                      threshHigh=0.8,doRel=TRUE,
                      doCat=FALSE,doScale=FALSE)}
\arguments{
\item{dir}{Directory with seg files in.}
\item{threshLow}{Threshold above which to call copy number gain or loss.}
\item{threshHigh}{Threshold above which to call extreme copy numebr gain or loss.}
\item{doRel}{Whether to adjust scores based on relative lengths of chromosomes.}
\item{doCat}{Whether to categorise copy number chaanges into chromosomal, arm-level, focal or large split.}
\item{doScale}{Whether to scale copy number scores.}
}
\value{
Matrix of copy number categories if \code{doCat=TRUE}, or vector of copy number scores if \code{doCat=FALSE}.
} 
\description{
Score chromosomal level copy number changes. 
}
\details{
  Implements copy number scores from Davoli et al. (2017).
}

