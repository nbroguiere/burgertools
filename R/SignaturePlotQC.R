#' QC plot with signatures highlight
#'
#' This function makes quality control (QC) plots from a Seurat object, by default percent.mito vs nFeature_RNA, while highlighting signature scores in color.
#' 
#' Can also be used in general as a combination of scatter plot and feature plot, as any metadata/signature/feature can be used for x and y axis as well as for color.
#'
#' This is useful when making scRNA-seq from heterogeneous samples: for example immune cells are small and might be discarded as debris from much larger epithelial cancer cells if one is not aware of the differences in QC parameters occurring between various cell types.
#' @param object Seurat object. Must contain nFeature_RNA and percent.mito metadata columns, or the x and y plot parameters must be changed accordingly to plot other QC info otherwise.
#' @param signatures character(n) or list(character(1)) or named list(n) vector containing the names of the signatures to highlight in color. Can be metadata columns or assay features. If a named list of signatures is passed, the names of the list will be used. If values are not found in the metadata, they will be picked from assay features. If the feature is found in several assays, it should be further defined with an underscore separator, in the form "assay_feature".
#' @param x Which parameter to set to the x axis in the QC plot (Default: nFeature_RNA).
#' @param y Which parameter to set to the y axis in the QC plot (Default: percent.mito).
#' @param log.scale logical(1). Whether to plot with log.scale or not (Default: TRUE)
#' @param ncol num(1). Number of columns (Default = NA, determined from the data)
#' @param pt.size num(1). The point size, passed to ggplot (Default: 1)
#' @param assay character(1). If some values to highlight are features rather than metadata columns, the assay from which they should be pulled in priority (Default: DefaultAssay(object)).
#' @param layer character(1). If some values to highlight are features rather than metadata columns, the layer from which they should be pulled (Default: "data").
#' @return Returns a ggplot/cowplot grid object.
#' @keywords QC plot signatures highlight
#' @export
#' @examples
#' # After MySeuratObject has been log-normalized in order to contain an RNA>data assay>layer and been augmented with a percent.mito metadata column:
#' # SignatureList should be a list of character vectors, each containing a series of feature/gene names.
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject <- ScoreSignatures(MySeuratObject,SignatureList)
#' SignaturePlotQC(MySeuratObject,names(SignatureList))
#' SignaturePlotQC(MySeuratObject,x="MKI67",y="PCNA",log.scale=FALSE,pt.size=1.5,)

SignaturePlotQC <- function(object, signatures, x= "nFeature_RNA", y="percent.mito", log.scale=TRUE, ncol=NA, pt.size=1, assay=DefaultAssay(object), layer="data"){
  if(is.list(signatures)){
    if(length(unlist(signatures))==length(signatures)){
      sign.names <- unlist(signatures)
    }else{
      sign.names <- names(signatures)
    }
  }else if(is.vector(signatures, mode="character")){
    sign.names <- signatures
  }else{
    stop("The format of the 'signatures' argument (i.e. ",typeof(signatures),") is not supported. Provide a vector of metadata column names, feature names, assay_feature names, or a named signature list.")
  }
  # sign.not.found <- setdiff(sign.names, colnames(object@meta.data))
  # if(length(sign.not.found)>0){
  #   cat(paste("Signatures not found:",toString(sign.not.found),"\nAttempt to use these features from imputed or default assay instead.\n"))
  #   sign.names <- intersect(sign.names, colnames(object@meta.data))
  # }
  # df <- object@meta.data
  # imputed_use <- intersect(rownames(object[["imputed.RNA"]]),sign.not.found)
  # if(length(imputed_use)){
  #   for(i in imputed_use)
  #     df[,i] <- as.numeric(object[["imputed.RNA"]][i,])
  # }
  # default_use <- intersect(rownames(object),setdiff(sign.not.found,imputed_use))
  # if(length(default_use)){
  #   for(i in default_use)
  #     df[,i] <- as.numeric(object[[DefaultAssay(object)]][i,])
  #   warning("The format of the 'signatures' argument (i.e. ",typeof(signatures),") is not supported. Provide a vector of metadata/feature/signature names, or a named signature list.")
  #   return()
  # }
  df <- GatherFeatures(object,c(x,y,sign.names), layer=layer, assay=assay)
  p <- list()
  for(n in sign.names){
    if(log.scale){
      p[[n]] <- ggplot(df) + geom_point(aes_string(x=x, y=y, color=n), size=pt.size) + scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10') + lims(colour=c(0,NA))
    }else{
      p[[n]] <- ggplot(df) + geom_point(aes_string(x=x, y=y, color=n), size=pt.size) + lims(colour=c(0,NA))
    }
  }
  if(is.na(ncol)){
    if(length(sign.names)==1){
      ncol <- 1
    }else if(length(sign.names)>9){
      ncol <- 4
    }else if(length(sign.names)>4){
      ncol <- 3
    }else if(length(sign.names)>1){
      ncol <- 2
    }else{
      warning(paste("Incorrect number of signatures found (",length(sign.names),")."))
    }
  }
  p <- cowplot::plot_grid(plotlist = p, ncol= ncol)
  return(p)
}

