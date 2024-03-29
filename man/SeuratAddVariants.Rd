% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeuratAddVariants.R
\name{SeuratAddVariants}
\alias{SeuratAddVariants}
\title{Add the variant information contained in a genotype object to a Seurat object}
\usage{
SeuratAddVariants(
  seurat,
  genotype,
  consensus = "VAR",
  reference = "REF",
  alternate = "ALT",
  frequency = "FREQ"
)
}
\arguments{
\item{seurat}{Seurat object.}

\item{genotype}{Genotype object.}

\item{consensus}{character(1). The name of the Seurat assay that will store the consensus matrix (values 0 = no data, 1 = ref/ref, 2=alt/alt, 3=ref/alt).}

\item{reference}{character(1). The name of the Seurat assay that will store the reference allele count matrix.}

\item{alternate}{character(1). The name of the Seurat assay that will store the alternate allele count matrix.}

\item{frequency}{character(1). The name of the Seurat assay that will store the frequency matrix (0 = no data, 1-2 = ref only to alt only).}
}
\description{
The genotype object should contain a consensus call (vartrix conventions), reference allele count, alternate allele count, and a frequency matrix (typically offset by one to distinguish 0 frequency from no data in an efficient sparse matrix format), as well as metadata about the variants. By default, the single cell matrices will be stored in assays named VAR, REF, ALT and FREQ respectively, and the metadata will be added to the VAR assay feature metadata.
}
\examples{
# MySeuratObject <- SeuratAddVariants(MySeuratObject, MyGenotypeObject)
}
\keyword{matrices}
\keyword{seurat}
\keyword{variant}
\keyword{vartrix}
