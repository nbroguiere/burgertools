#' Automatic cell type classification based on signatures
#'
#' This function performs a simple automatic cell type classification based on signatures or other scores stored in the metadata or features (unique name or with in the format assay_feature). Useful for automating the simple first broad classification of cell types before a cell-type aware QC (for example, taking into account that immune cells have less detected genes/transcripts, and fibroblast less mitochondrial content, than epithelial/cancer cells).
#'
#' Each single cell will be associated to the cell type whose signature is maximal on the cell. A list of expected values for the signatures is optional (typically obtained by checking the typical signature score for cells of known identity beforehand). If given, expected values are used to normalize signature scores before classifying cells to the signature having maximal score. It is strongly recommended to use imputed data to score cell type signatures, to avoid misclassifications due to dropouts.
#'
#' @param object Seurat object. Must contain nFeature_RNA and percent.mito metadata columns, or the x and y plot parameters must be changed accordingly to plot other QC info otherwise.
#' @param signatures character(n) or named list(n). Names of the signatures to use. Must correspond to names of metadata columns or features is the current default assay. If a named list of signatures is passed, the names of the list will be used and the content ignored.
#' @param expected.values numeric(n). Expected values for the signatures, used for normalization of the signatures before classifying cells. If NA (Default), normalized to max. If 1, no normalization. If The vector/list is named and the names match the signature names, the names will be used to match each scale factor to the right signature. If no matching names are given, the expected values are assumed to be given in the same order as the signatures.
#' @param metadata.name character(1). The name of the new metadata column where cell type annotations are stored (Default: celltype)
#' @param cell.names character(n). Give only if the signature names should not be used as celltypes names, but rather be replaced by these cell names.
#' @param assay character(1). If some signatures used for classification are stored as assay features, which assay should be used in priority (Default: DefaultAssay(object)).
#' @param slot character(1). If some signatures used for classification are stored as assay features, which slot should be used (Default: "data").
#' @param restrict.ident character(1). The name of a metadata column containing celltype annotations, which will be used to define the subset of cells that should be classified. Default: current Idents(object). 
#' @param restrict.to character(n). Which celltypes (as defined in the restrict.ident metadata column) should be classified. Default: All cells. 
#' @return A Seurat object with an additional metadata column containing the cell type annotations and Idents() set to these annotations.
#' @keywords Cell type classification celltype Classifier
#' @export
#' @examples
#' # After MySeuratObject has been log-normalized in order to contain an RNA>data assay>slot and been augmented with a percent.mito metadata column:
#' # SignatureList should be a list of character vectors, each containing a series of feature/gene names. In this example, the signatures correspond to cell types, and the name of each signature is the name of the corresponding cell type (e.g. Tcell, Bcell, EpithelialCell, Fibroblasts, etc).
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject <- ScoreSignatures(MySeuratObject,SignatureList)
#' MySeuratObject <- Classify(MySeuratObject,SignatureList) # Automatic annotation based on cell type signatures, stored in metadata column "celltype".

Classify <- function(object, signatures, expected.values=NA, metadata.name="celltype", cell.names=NA, assay=DefaultAssay(object), slot="data", restrict.ident="Default", restrict.to=as.character(unique(Idents(object)))){

  # Check that the requested Idents to use are correct
  if(!length(restrict.ident)){
    warning("Ident to use for restriction not found. Aborting.")
    return(object)
  }else{
    if(!restrict.ident %in% c("Default",colnames(object@meta.data))){
      warning("Ident to use for restriction not found. Aborting.")
      return(object)
    }
  }
  
  # Define which cells should be classified
  if(!length(restrict.to)){
    warning("The restrict.to argument is empty, which means no cells should be classified.")
    return(object)
  }else{
    if(restrict.ident == "Default"){
      cells.use <- colnames(object)[Idents(object) %in% restrict.to]
    }else{
      cells.use <- colnames(object)[object@meta.data[,restrict.ident,drop=T] %in% restrict.to]
    }
    n.cells <- length(cells.use)
    if(!n.cells){
      warning("There are no cells in the requested restriction. Check the celltypes required as restrict.to exist within the restrict.ident.")
      return(object)
    }
  }
  
  # If the metadata column does not exist, create it, filled with "undefined"
  if(!metadata.name %in% colnames(object@meta.data)){
    object@meta.data[,metadata.name] <- "undefined"
  }
  
  # Check the signatures
  if(is.list(signatures)){
    sign.names <- names(signatures)
  }else if(is.vector(signatures, mode="character")){
    sign.names <- signatures
  }else{
    warning("The format of the 'signatures' argument (i.e. ",typeof(signatures),") is not supported. Provide a vector of metadata column / current assay feature names, or a named signature list.")
    return(object)
  }

  # Create a df with the signatures and features to be used:
  df <- GatherFeatures(object, sign.names, assay=assay, slot=slot)
  df <- df[cells.use,,drop=F]
  
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

  # Normalize the signatures to their expected values.
  sign.scores <- data.frame(row.names = rownames(df))
  for(n in sign.names){
    sign.scores[,n] <- df[,n]/expected.values[n]
  }

  # Classify to max signature
  object@meta.data[cells.use,metadata.name] <- colnames(sign.scores)[max.col(m = sign.scores)]

  # Rename the cells if cell.names is defined:
  if(length(cell.names)){
    if(!is.na(cell.names[1])){
      cell.names <- setNames(cell.names,sign.names)
      object@meta.data[cells.use,metadata.name] <- cell.names[object@meta.data[cells.use,metadata.name]]
    }
  }
  
  # Make it the default Idents and return
  Idents(object) <- object[[metadata.name]]
  return(object)
}
