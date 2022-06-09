#' Prefilter variants by coverage
#'
#' Given a list of vartrix matrices, and a minimal number of cells/barcodes for which variant data should be available, returns a list of vartrix matrices with filtered variants. If barcodes are given, additionally restrict the matrices to these barcodes.
#' The row order and rownames are kept consistent with the input and across all matrices during the subsetting.
#'
#' @param vartrix list of dgCmatrix, a list of vartrix matrices (containg at least a consensus matrix, named as such, in order to compute coverage).
#' @param min.cells numeric(1) the minimum number of cells for which there should be a genotype call for the variants to be kept.
#' @param barcodes_keep character(n) A vector of cell barcodes to keep. NA to keep all. Default:NA.
#' @return A list of vartrix matrices in the same format as the input
#' @keywords genotypes variants vcf
#' @examples
#' MyVartrixMatrices <- PrefilterVartrix(MyVartrixMatrices, min.cells=10)
#' @export

PrefilterVartrix <- function(vartrix_matrices,min.cells=5,barcodes_keep=NA){
  if(!is.na(barcodes_keep[1])){
    vartrix_matrices <- lapply(vartrix_matrices,function(m) return(m[,barcodes_keep]))
  }

  coverage <- Matrix::rowSums(vartrix_matrices$consensus>0)
  keep <- coverage>=min.cells
  return(lapply(vartrix_matrices,"[",keep,T))
}
