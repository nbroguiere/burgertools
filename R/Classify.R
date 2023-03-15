#' Automatic cell type classification based on signatures
#'
#' This function performs a simple automatic cell type classification based on signatures or other scores stored in the metadata or features (features names should be unique accross assays, or in the format assay_feature). Useful for automating the simple first broad classification of cell types before a cell-type aware QC (for example, taking into account that immune cells have less detected genes/transcripts, and fibroblast less mitochondrial content, than epithelial/cancer cells). Also useful in general for classifying cells to their most likely identity, while not attempting to remove any multiplet or cell with no good identity match (for example to look for subcategories such as polarities within a previously identified celltype). 
#'
#' Each single cell will be associated to the cell type whose signature is maximal on the cell. A list of expected values for the signatures is optional (typically obtained by checking the typical signature score for cells of known identity beforehand). If given, expected values are used to normalize signature scores before classifying cells to the signature having maximal score. It is strongly recommended to use imputed data to score cell type signatures, to avoid misclassifications due to dropouts.
#'
#' If the by.clusters parameter is set to some identities or metadata column containing identities, the classification will be done on clusters rather than single cells. This can help to overcome strong noise, or to take into account some prior knowledge such as dataset alignment in the classification.
#'
#' @param object Seurat object. Must contain nFeature_RNA and percent.mito metadata columns, or the x and y plot parameters must be changed accordingly to plot other QC info otherwise.
#' @param signatures character(n) or named list(n). Names of the signatures to use. Must correspond to names of metadata columns or features is the current default assay. If a named list of signatures is passed, the names of the list will be used and the content ignored.
#' @param expected.values numeric(n). Expected values for the signatures, used for normalization of the signatures before classifying cells. If NA (Default), normalized to max. If 1, no normalization. If The vector/list is named and the names match the signature names, the names will be used to match each scale factor to the right signature. If no matching names are given, the expected values are assumed to be given in the same order as the signatures.
#' @param metadata.name character(1). The name of the new metadata column where cell type annotations are stored (Default: celltype)
#' @param cell.names character(n). Give only if the signature names should not be used as celltypes names, but rather be replaced by these cell names.
#' @param assay character(1). If some signatures used for classification are stored as assay features, which assay should be used in priority (Default: DefaultAssay(object)).
#' @param slot character(1). If some signatures used for classification are stored as assay features, which slot should be used (Default: "data").
#' @param restrict.ident character(1). The name of a metadata column containing celltype annotations, which will be used to define the subset of cells that should be classified. Default: "Default", which uses current default, i.e. Idents(object). 
#' @param restrict.to character(n). Which celltypes (as defined in the restrict.ident metadata column) should be classified. Default: All cells. 
#' @param by.clusters character(1) or character(n) or factor(n). If not NULL, clusters instead of single cells will be classified. The parameter should indicate the name of a metadata column containing the clusters to consider (e.g. "seurat_clusters"), or the actual cluster info (e.g. object$seurat_clusters). 
#' @return A Seurat object with an additional metadata column containing the cell type annotations and Idents() set to these annotations.
#' @keywords Cell type classification celltype Classifier
#' @export
#' @examples
#' # After MySeuratObject has been log-normalized in order to contain an RNA>data assay>slot and been augmented with a percent.mito metadata column:
#' # SignatureList should be a list of character vectors, each containing a series of feature/gene names. In this example, the signatures correspond to cell types, and the name of each signature is the name of the corresponding cell type (e.g. Tcell, Bcell, EpithelialCell, Fibroblasts, etc).
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject <- ScoreSignatures(MySeuratObject,SignatureList)
#' MySeuratObject <- Classify(MySeuratObject,SignatureList) # Automatic annotation based on cell type signatures, stored in metadata column "celltype".

Classify <- function(object, signatures, expected.values=NA, metadata.name="celltype", cell.names=NA, assay=DefaultAssay(object), slot="data", restrict.ident="Default", restrict.to=as.character(unique(Idents(object))), by.clusters=NULL){
  
  # Check that the idents requested for restriction are correct
  if(length(restrict.ident) != 1){
    stop("restrict.ident should be of length 1. No classification was done.")
  }else{
    if(!restrict.ident %in% c("Default",colnames(object@meta.data))){
      stop("Ident to use for restriction not found. No classification was done.")
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
  
  # Check whether or not the classification should be done by clusters, and if yes determine what type of argument was given.
  if(length(by.clusters)){
    if(length(by.clusters)==1){ # Name of a metadata column was given.
      metadata.column.name <- by.clusters
      clusters <- as.character(object@meta.data[cells.use,metadata.column.name])
      print(paste0("Classifying by cluster, using the metadata column: ", metadata.column.name))
    }else if(length(by.clusters)==ncol(object)){ # complete list of clusters was given.
      clusters <- as.character(by.clusters)[colnames(object) %in% cells.use]
      print("Classifying by cluster.")
    }else{
      stop("The length of the by.clusters argument should be 1 (metadata column name) or n cells (vector of idents for all cells).")
    }
  }else{
    print("No clusters given, classifying single cells.")
  }
  
  # If the metadata column for the output does not exist, create it, filled with "undefined"
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
  if(!ncol(df)){
    warning("None of the requested features/signatures were found. Did you forget to run Impute and/or ScoreSignatures?")
    return(object)
  }
  
  # If classifying by clusters, average out df over clusters:
  if(length(by.clusters)){
    unique.clusters <- unique(clusters)
    df.avg <- as.data.frame(matrix(0, nrow = length(unique.clusters), ncol = length(sign.names), dimnames = list(unique.clusters,sign.names)))
    for(i in unique.clusters){
      cells.tmp <- rownames(df)[clusters==i]
      for(j in sign.names){
        df.avg[i,j] <- mean(df[cells.tmp,j])
      }
    }
    df <- df.avg
  }
  
  # Determine which mode should be used for expected.values
  if(sum(is.na(expected.values))){
    print("expected.values=NA - Normalizing the signatures to their max.")
    expected.values = MatrixGenerics::colMaxs(as.matrix(df))
    names(expected.values) <- sign.names
    print(expected.values)
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
  sign.scores <- t(t(df)/expected.values[colnames(df)])


  # Classify to max signature
  if(!length(by.clusters)){ # By single cells
    object@meta.data[cells.use,metadata.name] <- colnames(sign.scores)[max.col(m = sign.scores)]
  }else{ # By cluster
    assignments <- setNames(colnames(df)[apply(sign.scores, 1, which.max)], rownames(sign.scores))
    object@meta.data[cells.use,metadata.name] <- assignments[clusters]
  }
  
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
