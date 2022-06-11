#' Prefilter variants by coverage in the vartrix matrices
#'
#' Given a genotype object containing a list of vartrix matrices, and a minimal number of cells/barcodes for which variant data should be available, returns a list of vartrix matrices with filtered variants. If barcodes are given, start by restricting the matrices to these barcodes before filtering variants.
#' The row order and rownames are kept consistent with the input during the subsetting. The subsetting is also applied to the entire genotype object, including reference genotypes and variant metadata. Stores the coverage information in genotype$coverage.
#'
#' @param genotype A genotype object, containing vartrix matrices (including a VAR consensus matrix, named as such, used to compute coverage).
#' @param min.cells numeric(1) the minimum number of cells for which there should be a genotype call for the variants to be kept.
#' @param barcodes_keep character(n) A vector of cell barcodes to keep. NA to keep all. Default:NA.
#' @return The filtered genotype object
#' @keywords genotypes variants vcf
#' @examples
#' MyFilteredGenotypeObject <- PrefilterVartrix(MyGenotypeObject, min.cells=10)
#' @export

PrefilterVartrix <- function(genotype, min.cells=5, barcodes_keep=NA){
  # Subset the barcodes first if applicable:
  if(!is.na(barcodes_keep[1])){
    genotype@vartrix <- lapply(genotype@vartrix,function(m) return(m[,barcodes_keep]))
  }

  genotype$coverage <- Matrix::rowSums(genotype@vartrix$VAR>0)
  keep <- genotype@variants[genotype$coverage>=min.cells]
  genotype@vartrix <- lapply(genotype@vartrix,"[",keep,T)
  return(genotype)
}
