#' Compute the fraction/percentage of mitochondrial genes (mouse and human) in Seurat objects
#'
#' This function computes the mitochondrial content of cells in a Seurat object and adds it to the metadata.
#' @param object Seurat object. Must contain an RNA assay.
#' @param name Name of the new metadata column created containing the fraction of mitochondrial transcripts (Default: percent.mito)
#' @param as.percent Whether to report the result as a fraction (between 0 and 1) or percentage (between 0 and 100). Default: TRUE (fraction).
#' @return Seurat object with an additional metadata column containing mitochondrial fraction/percentage.
#' @keywords Mitochondrial genes content percentage fraction Mito
#' @export
#' @examples
#' MySeuratObject <- ComputeMitoContent(MySeuratObject)

ComputeMitoContent <- function(object, name = "percent.mito", as.percent = FALSE){

  mouse_mito_genes <- grep(pattern = "^mt-", x = rownames(object), value = TRUE)
  human_mito_genes <- grep(pattern = "^MT-", x = rownames(object), value = TRUE)
  n_mm <- length(mouse_mito_genes)
  n_hu <- length(human_mito_genes)
  if(max(n_mm,n_hu)==0){
    warning("No mitochondrial genes (starting with mt- or MT-) were found in the Seurat object. No metadata column added.")
    return(object)
  }else{
    if(n_hu>n_mm){
      mito.content <- colSums(object[grep(pattern = "^MT-", x = rownames(object), value = TRUE), ])/colSums(object)
    }else{
      mito.content <- colSums(object[grep(pattern = "^mt-", x = rownames(object), value = TRUE), ])/colSums(object)
    }
    if(as.percent){
      mito.content <- 100*mito.content
    }
    object <- AddMetaData(object = object, metadata = mito.content, col.name = name)
    return(object)
  }
}
