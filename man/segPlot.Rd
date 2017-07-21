\name{segPlot} 
\alias{segPlot} 
\title{CNomplexity::segPlot} 
\usage{segPlot(segFile=NULL,segObj=NULL,dataDir,
		outDir=NULL,fileName=NULL,
		sampleCol=2,chromCol=1,startCol=3,endCol=4,nMarkCol=NULL,segMeanCol=6,
		plotAber=FALSE,highlightChroms=NULL,highlightStarts=NULL,highlightEnds=NULL,
		fusionChroms=NULL,fusionPos=NULL,segColour=NULL)}
\arguments{
\item{segFile}{Segment file to plot.}
\item{segObj}{Segment object to plot.}
\item{dataDir}{Data directory if \code{segFile!=NULL}.}
\item{outDir}{Out directory for plotting.}
\item{fileName}{OUt filename for plotting.}
\item{sampleCol}{Sample column of seg profile.}
\item{chromCol}{Chromsome column of seg profile.}
\item{startCol}{Start column of seg profile.}
\item{endCol}{End column of seg profile.}
\item{nMarkCol}{Seg column with number of markers.}
\item{segMeanCol}{Seg column with copy number.}
\item{highlightChroms}{Chromsomes of regions to highlight.}
\item{highlightStarts}{Start positions of regions to highlight.}
\item{highlightEnds}{End positions of regions to highlight.}
\item{fusionChroms}{Chromosomes of fusions to plot.}
\item{fusionPos}{Positions of fusions to plot.}
\item{segColour}{Colour for segments.}
}
\value{
NULL
} 
\description{
Function to plot a segment profile., with options to higlight some regions, and plot vertical bars at fusion positions.
}
\details{
  Plots a segment profile, with options to higlight some regions, and plot vertical bars at fusion positions.
}

