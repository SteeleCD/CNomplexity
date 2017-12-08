name{nmfFitMetrics} 
\alias{nmfFitMetrics} 
\title{CNomplexity::nmfFitMetrics} 
\usage{nmfFitMetrics(N,res,data,
	plotsil=FALSE,maxiter=1000,nstarts=20,
	algo="Hartigan-Wong",doMetric=TRUE,
	clusterMeth="kmeans",hclustInput="sigs",
	kmeansDist="cosine") }
\arguments{
\item{N}{String. Which rank of \code{res} to generate metrics for.}
\item{res}{Results object from a \code{replicate} of \code{bootNMF}.}
\item{data}{Matrix of n*m that was handed to bootNMF. n=number of samples, m=number of features.}
\item{plotsil}{Boolean. Whetehr to plot results of silhouette.}
\item{maxiter}{Maximum iterations for kmeans. Only used if \code{kmeansDist} is not \code{cosine}.}
\item{nstarts}{Number of starts for kmeans. Only used if \code{kmeansDist} is not \code{cosine}.}
\item{algo}{Algorithm for kmeans. Only used if \code{kmeansDist} is not \code{cosine}.}
\item{doMetric}{Boolean. Whether to calculate metrics.}
\item{clusterMeth}{Cluster method. Defaults to k-means clustering. Otherwise, heirarchical clustering.}
\item{hclustInput}{Input to heirarchical clustering. Default is the raw signatures. Otherwise, cosine distances between signatures.}
\item{kmeansDist}{Distance metric for kmeans. Defaults to cosine. Otherwise, euclidean distance.}
}
\value{
If \code{doMetric}=\code{TRUE}, a list with 4 elements. \code{sil}=average silhouette width across groups. \code{frob}=Frobenius distance between the reconstructed data from signatures and the original data. \code{P}=extracted signatures, m*N matrix where m is the number of features, and N is the number of signatures to extract. \code{E}=exposures, N*n matrix, where n is the numebr of samples in original data. If \code{doMetric}=\code{FALSE}, only \code{P} and \code{E} are returned. 
} 
\description{
Cluster signatures extracted from multiple runs of nmfBoot for a given number of signatures, N, into N clusters. Calculate Frobenius distance between reconstructed data from the means of these clustered signatures and the original data, as well as the average silhouette width of the members of each clustered signature. Aim to maximise the silhouette width, while minimising the Frobenius distance. 
}
\details{
  Method for selecting the number of signatures that are appropriate for a given dataset. Described in Alexandrov et al. (2013).
}

