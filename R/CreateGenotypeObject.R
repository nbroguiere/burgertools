#' Create a genotype object
#'
#' @param matrix A dgCMatrix containing genotype data in vartrix conventions.
#' @param metadata A data.frame containing variant metadata.
#' @param variants A character vector listing the variants described in the object, in the same order as the rows of the matrix and metadata slots.
#' @param vartrix A list of vartrix matrices
#' @param variants_by_coverage A character vector listing the variants sorted from max coverage (number of cells in which there is data) to min.
#' @param variants_by_information A character vector listing the variants sorted from max information (excess entropy in single cell data) to min.
#' @param informative_variants A character vector listing the most informative variants, used for downstream clustering analysis.
#' @return Returns a genotype object.
#' @keywords genotype vcf vep variants
#' @export
CreateGenotypeObject <- function(matrix=NULL, metadata=NULL, vartrix=NULL, variants=NULL){
  if(is.null(variants)){
    if(rownames(matrix)!=rownames(metadata)){
      warning("The row names in the matrix and metadata do not match, and list of variants is not given. Aborting.")
    }else{
      variants <- rownames(matrix)
    }
  }
  rownames(matrix) <- variants
  rownames(metadata) <- variants
  rownames(vartrix) <- variants
  return(GenotypeObject(matrix=matrix,metadata=metadata,vartrix=vartrix,variants=variants))
}
