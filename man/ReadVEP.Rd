% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadVep.R
\name{ReadVEP}
\alias{ReadVEP}
\title{Read a VEP tsv file}
\usage{
ReadVEP(file, remove_empty_col = TRUE)
}
\arguments{
\item{file}{character(1) The path to the VEP file, can be compressed (.vcf or .vcf.gz).}

\item{remove_empty_col}{logical(1) Should uninformative columns with only a unique value/level be removed. Default: TRUE.}
}
\value{
A tibble
}
\description{
This function reads a variant effect predictor (VEP) output file in tsv format, and returns the contents as a tibble.
}
\details{
This function is a readr::read_tsv wrapper that reads a variant effect predictor (VEP) output file in tsv format, and returns the contents as a tibble. By default, removes uninformative columns (i.e. which have only one value), and setups the IMPACT (if present) as an ordered factor.
}
\examples{
genotypes <- ReadVEP("MyFolder/MyVEPresults.tsv.gz")
}
\keyword{genotypes}
\keyword{variants}
\keyword{vcf}
