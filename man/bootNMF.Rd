\name{bootNMF} 
\alias{bootNMF} 
\title{CNomplexity::bootNMF} 
\usage{bootNMF(data,rank) }
\arguments{
\item{data}{Matrix of n*m, where n is the number of samples, and m is the number features.}
\item{rank}{Integer vector. Vector of the number of signatures to identify.}
}
\value{
List wih two elements. \code{nmf} = nmf results for each rank. \code{sil} = silhouette results for each rank.
} 
\description{
Bootstrap input data, and run non-negative matrix factorisation on the bootstrapped data to identify latent signatures in the original data. Should be run multiple times e.g. in a \code{replicate} function.
}
\details{
  Method for bootstrapping and running non-negative matrix factorisation on a dataset, as described in Alexandrov et al. (2013).
}

