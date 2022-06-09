#' Automatic cell type classification based on signatures
#'
#' This function performs a simple automatic cell type classification based on signatures or other scores stored in the metadata or features in the current active assay, data slot. Useful for automating the simple first broad classification of cell types before a cell-type aware QC (for example, taking into account that immune cells have less detected genes/transcripts, and fibroblast less mitochondrial content, than epithelial/cancer cells).
#'
#' Each single cell will be associated to the cell type whose signature is maximal on the cell. A list of expected values for the signatures is optional (typically obtained by checking the typical signature score for cells of known identity beforehand). If given, expected values are used to normalize signature scores before classifying cells to the signature having maximal score. It is strongly recommended to use imputed data to score cell type signatures, to avoid misclassifications due to dropouts.
#'
#' @param object Seurat object. Must contain nFeature_RNA and percent.mito metadata columns, or the x and y plot parameters must be changed accordingly to plot other QC info otherwise.
#' @param signatures character(n) or named list(n). Names of the signatures to use. Must correspond to names of metadata columns or features is the current default assay. If a named list of signatures is passed, the names of the list will be used.
#' @param expected.values numeric(n). Expected values for the signatures, used for normalization of the signatures before classifying cells. If NA (Default), normalized to max. If 1, no normalization. If The vector/list is named and the names match the signature names, the names will be used to match each scale factor to the right signature. If no matching names are given, the expected values are assumed to be given in the same order as the signatures.
#' @param metadata.name character(1). The name of the new metadata column where cell type annotations are stored (Default: celltype)
#' @return A Seurat object with an additional metadata column containing the cell type annotations and Idents() set to these annotations.
#' @keywords Cell type classification celltype Classifier
#' @export
#' @examples
#' # After MySeuratObject has been log-normalized in order to contain an RNA>data assay>slot and been augmented with a percent.mito metadata column:
#' # SignatureList should be a list of character vectors, each containing a series of feature/gene names. In this example, the signatures correspond to cell types, and the name of each signature is the name of the corresponding cell type (e.g. Tcell, Bcell, EpithelialCell, Fibroblasts, etc).
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject <- ScoreSignatures(MySeuratObject,SignatureList)
#' MySeuratObject <- Classify(MySeuratObject,SignatureList) # Automatic annotation based on cell type signatures, stored in metadata column "celltype".

Classify <- function(object, signatures, expected.values=NA, metadata.name="celltype"){
  if(is.list(signatures)){
    sign.names <- names(signatures)
  }else if(is.vector(signatures, mode="character")){
    sign.names <- signatures
  }else{
    warning("The format of the 'signatures' argument (i.e. ",typeof(signatures),") is not supported. Provide a vector of metadata column / current assay feature names, or a named signature list.")
    return(object)
  }
  sign.available <- c(colnames(object@meta.data),rownames(object))
  sign.not.found <- setdiff(sign.names, sign.available)
  if(length(sign.not.found)>0){
    warning(paste("Signatures/features not found, proceeding without them:",toString(sign.not.found)))
    sign.names <- intersect(sign.names, colnames(object@meta.data))
  }

  # Create a df with the signatures and features to be used:
  df <- cbind(object@meta.data[,intersect(sign.names,colnames(object@meta.data))],
              t(as.data.frame(object[[DefaultAssay(object)]]@data[setdiff(sign.names,colnames(object@meta.data)),])))

  if(sum(is.na(expected.values))){
    print("expected.values=NA - Normalizing the signatures to their max.")
    expected.values = MatrixGenerics::colMaxs(as.matrix(df))
    names(expected.values) <- sign.names
  }else if(length(expected.values)==1 & expected.values[1]==1){
    expected.values <- vector(mode="numeric",length=length(sign.names))+1
    names(expected.values) <- sign.names
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

  # Initialize the cell type filter (start from keep all, to then restrict iteratively for each condition):
  filters <- matrix(data = TRUE, nrow = dim(object)[2], ncol = length(celltype_names), dimnames = list(rownames(df),celltype_names))

  sign.scores <- data.frame(row.names = rownames(df))
  for(n in sign.names){
    sign.scores[,n] <- df[,n]/expected.values[n]
  }
  object[[metadata.name]] <- colnames(sign.scores)[max.col(m = sign.scores)]
  Idents(object) <- object[[metadata.name]]
  return(object)
}
