#' Create a genotype object
#'
#' The genotype objects are tailored to handle efficiently in R the contents of (mouse/human) vcf files, i.e. variant metadata, and several genotypes. They are also meant to be augmented with additional information such as annotations (e.g. vep) and single-cell analysis statistics.
#'
#' For efficient handling, the genotypes are stored as a sparse numeric matrix. When creating the object with ReadVcf, the vartrix conventions (from 10X genomics) are used, namely:
#' 0 for "no call" or in vcf "./."
#' 1 for "ref/ref" or in vcf "0/0"
#' 2 for "alt/alt" or in vcf "1/1"
#' 3 for "ref/alt" or in vcf "0/1"
#'
#' The objects can instead be created manually with CreateGenotypeObject().
#'
#' The genotype data is found in the matrix slot, and information concerning the variants is contained in the metadata slot. The list of variants currently in the object is found in the variants slot.
#'
#' Additional slots can be populated with additional data relevant to single cell analysis: most informative variants (informative_variants slot), and variants sorted by coverage (variants_by_coverage slot) or by entropy (variants_by_information).
#'
#' Standard generic methods can be used on the genotype object. In particular, the object can be subset with object[[i,j]] to obtain a genotype object with the corresponding subset of variants and genotypes (with full associated metadata, and restricted ranked variant lists). Object[i,j(,drop)] enables to access the genotype matrix in read and write, as well as rowSums, colSums, rownames, colnames, nrow, ncol, and dim. The $ operator enables to access directly to metadata columns, for both read and write.
#'
#' @param matrix A dgCMatrix containing genotype data in vartrix conventions.
#' @param metadata A data.frame containing variant metadata.
#' @param vartrix A list of vartrix matrices
#' @param variants A character vector listing the variants described in the object, in the same order as the rows of the matrix and metadata slots. Leave NULL to use the rownames of the matrix and metadata (They should match though).
#' @param variants_by_coverage A character vector listing the variants sorted from max coverage (number of cells in which there is data) to min.
#' @param variants_by_information A character vector listing the variants sorted from max information (excess entropy in single cell data) to min.
#' @param informative_variants A character vector listing the most informative variants, used for downstream clustering analysis.
#' @return Returns a genotype object.
#' @keywords genotype vcf vep variants
#' @export
CreateGenotypeObject <- function(matrix=NULL, metadata=NULL, vartrix=list(), variants=NULL, ...){
  if(is.null(variants)){
    if(!dplyr::setequal(rownames(matrix), rownames(metadata))){
      warning("The row names in the matrix and metadata do not match, and list of variants is not given. Aborting.")
    }else{
      variants <- rownames(matrix)
    }
  }
  if(!is.null(matrix)) rownames(matrix) <- variants
  if(!is.null(metadata)) rownames(metadata) <- variants
  if(length(vartrix)) rownames(vartrix) <- variants
  return(GenotypeObject(matrix=matrix, metadata=metadata, vartrix=vartrix, variants=variants))
}
