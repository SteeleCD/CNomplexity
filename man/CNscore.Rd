\name{CNscore} 
\alias{CNscore} 
\title{CNomplexity::CNscore} 
\usage{CNscore(segFiles,threshLow=0.275,
                      threshHigh=0.8,doRel=TRUE,
                      doCat=FALSE,doScale=FALSE,
			tumTotCol=7,normTotCol=5,
			startCol=3,endCol=4,chromCol=2,
			sampleCol=1,head=TRUE)}
\arguments{
\item{segFiles}{String or vector of strings. Path to segment files.}
\item{threshLow}{Numeric. Threshold above which to call copy number gain or loss.}
\item{threshHigh}{Numeric. Threshold above which to call extreme copy numebr gain or loss.}
\item{doRel}{Boolean. Whether to adjust scores based on relative lengths of chromosomes.}
\item{doCat}{Boolean. Whether to categorise copy number chaanges into chromosomal, arm-level, focal or large split.}
\item{doScale}{Boolean. Whether to scale copy number scores.}
\item{tumTotCol}{Integer. Seg column with total copy number information for tumour sample.}
\item{normTotCol}{Integer. Seg column with total copy number information for matched normal sample.}
\item{startCol}{Integer. Seg column with segment start positions.}
\item{endCol}{Integer. Seg column with segment end positions.}
\item{chromCol}{Integer. Seg column with segment chromosomes.}
\item{sampleCol}{Integer. Seg column with sample names for each segment.}
\item{head}{Boolean. Whether seg file has headers.}
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

