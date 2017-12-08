name{getRearrClass} 
\alias{getRearrClass} 
\title{CNomplexity::getRearrClass} 
\usage{getClustered(sample,chrom1,start1,end1,
	chrom2,start2,end2,type) }
\arguments{
\item{sample}{String vector. Sample for each rearrangement.}
\item{chrom1}{String vector. First partner chromosome for each rearrangement.}
\item{start1}{String vector. First partner start position for each rearrangement.}
\item{end1}{String vector. First partner end position for each rearrangement.}
\item{chrom2}{String vector. Second partner chromosome for each rearrangement.}
\item{start2}{String vector. Second partner start position for each rearrangement.}
\item{end2}{String vector. Second partner end position for each rearrangement.}
\item{type}{String vector. Structural variant type for each rearrangement e.g. deletion, inversion, tandem-duplication or translocation.}
}
\value{
Vector of rearrangement classes.
} 
\description{
Categorise rearrangements into classes for input into non-negative matrix factorisation. Categorised by rearrangement type, rearrangement size, and un/clu
stered status.
}
\details{
  Rearrangement classes are described in Nik-Zainal et al. (2016).
}

