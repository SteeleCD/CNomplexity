\name{chromothripsis} 
\alias{chromothripsis} 
\title{CNomplexity::chromothripsis} 
\usage{chromothripsis(segFile, bedpeFile=NULL, size=1e7,
			gap=1e6,chromCol=2,startCol=3,endCol=4,
			bedpeChromCol1=1,bedpePosCol1=2,direction1col=9,
			bedpeChromCol2=4,bedpePosCol2=5,direction2col=10,
			segSampleCol=1,bedpeSampleCol=1,doParallel=FALSE,
			nCores = NULL,samplesToRun = NULL,
			chromsToRun = NULL,sepbedpe=TRUE,
			bedpeHead=FALSE,segHead=TRUE,bedpeEnding=".brass.annot.bedpe.gz",
			breaksLimit=30,pThresh=0.8)}
\arguments{
\item{segFile}{Segment file.}
\item{bedpeFile}{Directory of separate bedpe files, or filename of single bedpe file.}
\item{size}{Window size for sliding window across chomosomes.}
\item{gap}{Gap between consecutive sliding windows.}
\item{chromCol}{Column of seg file with chromosome name.}
\item{startCol}{Column of seg file with start position.}
\item{endCol}{Column of seg file with end position.}
\item{bedpeChromCol1}{Column of bedpe file with chromosome name for first fusion partner.}
\item{bedpePosCol1}{Column of bedpe file with position for first fusion partner.}
\item{direction1col}{Column of bedpe file with orientation for first fusion partner.}
\item{bedpeChromCol2}{Column of bedpe file with chromosome name for second fusion partner.}
\item{bedpePosCol2}{Column of bedpe file with position for second fusion partner.}
\item{direction2col}{Column of bedpe file with orientation for second fusion partner.}
\item{segSampleCol}{Column of seg file with sample name.}
\item{bedpeSampleCol}{Column of bedpe file with sample name.}
\item{doParallel}{Whether to run in parallel}
\item{nCores}{Number of cores if \code{doParallel=TRUE}.}
\item{samplesToRun}{Character vector with sample names to run analysis for.}
\item{chromsToRun}{Character vector with chromosome names to run analysis for.}
\item{sepbedpe}{Whether bedpe files are separate or concatenated.}
\item{bedpeHead}{Whether bedpe file has a header.}
\item{segHead}{Whether seg file has a header}
\item{bedpeEnding}{File ending for separate bedpe files.}
\item{breaksLimit}{Minimum number of breakpoints on a chromsome to consider chromothripsis. Defaults to 30.}
\item{pThresh}{P value threshold for test of randomness (random order of breakpoints and random joins). Defaults to 0.8.}
}
\value{
List of samples and chromosmores, with each element being a matrix with Sample name, start positon of chromothriptic region and end position of chromothriptic region. NULL is no chromothriptic regions in the chromosome.
} 
\description{
Identify chromothriptic regions. 
}
\details{
  Implements three tests from Korbel & Campbell (2013).
}

