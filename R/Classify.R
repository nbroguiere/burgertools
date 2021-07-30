#' Cell classification based on signatures
#'
#' This function performs a simple cell type classification based on signatures, useful for automating the simple first broad classification of cell types before a cell-type aware QC (for example, taking into account that immune cells have less detected genes/transcripts, and fibroblast less mitochondrial content, than epithelial cancer cells).
#'
#' Each single cell will be associated to the cell type whose signature is maximal on the cell. A list of expected values for the signatures can be given (for example, checking the typical signature score for cells of known identity beforehand) and are used to normalize signature scores before determining the classification. It is strongly recommended to use imputed data to score the cell type signatures (markers), to avoid mis-classification due on dropouts.
#'
#' @param object Seurat object. Must contain nFeature_RNA and percent.mito metadata columns, or the x and y plot parameters must be changed accordingly to plot other QC info otherwise.
#' @param sign.names The names of the signatures to use for classification. Must correspond to the names of metadata columns. These names will also be used to name the corresponding cell types.
#' @param expected.values Vector/List of expected values for the signatures, used for normalization of the signatures before classifying cells. If NA (Default), no normalization is performed. If The vector/list is named and the names match the signature names, the names will be used to match each scale factor to the right signature. If no matching names are given, the expected values are assumed to be given in the same order as the signatures.
#' @param metadata.name The name of the new metadata column where cell type annotations are stored (Default: celltype)
#' @return A Seurat object with an additional metadata column containing the cell type annotations and Idents() set to these annotations.
#' @keywords Cell type classification celltype Classifier
#' @export
#' @examples
#' # After MySeuratObject has been log-normalized in order to contain an RNA>data assay>slot and been augmented with a percent.mito metadata column:
#' # SignatureList should be a list of character vectors, each containing a series of feature/gene names. In this example, the signatures correspond to cell types, and the name of each signature is the name of the corresponding cell type (e.g. Tcell, Bcell, EpithelialCell, Fibroblasts, etc).
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject <- ScoreSignatures(MySeuratObject,SignatureList)
#' MySeuratObject <- Classify(MySeuratObject,names(SignatureList),"CellType") # Automatic cell type annotation based on cell type signatures.

Classify <- function(object, sign.names, expected.values=NA, metadata.name="celltype"){
  if(sum(is.na(expected.values))){
    expected.values <- vector(mode="numeric",length=length(sign.names))+1
  }
  sign.not.found <- setdiff(sign.names, colnames(object@meta.data))
  if(length(sign.not.found)>0){
    warning(paste("Signatures not found:",toString(sign.not.found)))
    sign.names <- intersect(sign.names, colnames(object@meta.data))
  }
  if(length(expected.values)!=length(sign.names)){
    warning("The length of expected values does not match the length of the signature names. Aborting.")
    return(object)
  }
  if(length(intersect(names(expected.values),sign.names))!=length(sign.names)){
    if(length(intersect(names(expected.values),sign.names))==0){
      print("No matching names found in the expected values vector, assuming values are given in matching order.")
      names(expected.values) <- sign.names
    }else{
      warning("The names in the expected values vector do not match the ones in the signature names. Aborting.")
      return(object)
    }
  }
  sign.scores <- data.frame(row.names = rownames(object@meta.data))
  for(n in sign.names){
    sign.scores[,n] <- object@meta.data[,n]/expected.values[n]
  }
  object[[metadata.name]] <- colnames(sign.scores)[max.col(m = sign.scores)]
  Idents(object) <- object[[metadata.name]]
  return(object)
}
