% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindInformativeVariants.R
\name{FindInformativeVariants}
\alias{FindInformativeVariants}
\title{Find informative variants}
\usage{
FindInformativeVariants(
  seurat,
  genotype = NA,
  n.variants = 10000,
  assay = "VAR"
)
}
\arguments{
\item{seurat}{A seurat object containing a consensus matrix.}

\item{genotype}{A genotype object to be updated with the coverage and information data, or NA to give no genotype to update and return only a Seurat object.}

\item{n.variants}{numeric(1). The number of most informative variants stored in the @informative_variants slot (Default: 10000).}

\item{assay}{character(1). The name of the assay which stores the single cell variant calls data within the Seurat object, in vartrix formatting (sparse matrix with values 1-2-3 for ref-alt-het).}
}
\value{
Returns the updated Seurat object if no genotype is given, or list(seurat,genotype) objects if a genotype is given.
}
\description{
Computes the fraction of detection of variants at the single cell level, and entropy associated with the variants, based on a seurat objects containing a variant call consensus matrix (VAR). Stores this information in the assay metadata, and compute the most informative variants (stored as VariableFeatures of the assay). Returns the updated Seurat object. If a genotype object is also given, it will also be updated with this info, as well as variants_by_coverage and variants_by_information sorted lists.
}
\details{
FindInformativeVariants populates the following metadata columns and slots in the genotype/Seurat objects:

$excess_entropy (see below),

$coverage (number of cells with data for the variant),

$coverage_frac (fraction of the cell with data for the variant),

@variants_by_coverage (variants sorted from top to least coverage, genotype object only),

@variants_by_information (variants sorted from top to least excess entropy, genotype object only),

@informative_variants or VariableFeatures() (most informative variants, categorical data equivalent of what most variable features is for continuous data).

The amount of information provided by a variant at the single cell level is computed as the Shannon entropy, sum(-p*log2(p)), minus the minimal entropy contributed by the mere coverage of the variant, i.e. the entropy value if the data had just the two levels nodata and data. This is analogous to the excess variance for continuous data, translated to discrete variant data. This quantity is referred to as excess_entropy in the genotype metadata.
}
\examples{
MySeuratObject <- FindInformativeVariants(MySeuratObject, n.variants = 20000)
library(zeallot) # To enable the multiassignment operator. Otherwise need to deconstruct the list manually.
c(MySeuratObject, MyGenotypes) \%<-\% FindInformativeVariants(MySeuratObject, MyGenotypes, n.variants = 20000)
}
\keyword{entropy}
\keyword{excess}
\keyword{features}
\keyword{genotype}
\keyword{informative}
\keyword{most}
\keyword{variable}
\keyword{variants}
