#' QC plot with signatures highlight
#'
#' This function makes quality control (QC) plots from a Seurat object, by default percent.mito vs nFeature_RNA, while highlighting signature scores in color.
#'
#' This is useful when making scRNA-seq from heterogeneous samples: for example immune cells are small and might be discarded as debris from much larger epithelial cancer cells if one is not aware of the differences in QC parameters occuring between various cell types.
#' @param object Seurat object. Must contain nFeature_RNA and percent.mito metadata columns, or the x and y plot parameters must be changed accordingly to plot other QC info otherwise.
#' @param signatures A vector containing the names of the signatures to highlight in color. Must correspond to names of metadata columns. If a named list of signatures is passed, the names of the list will be used.
#' @param x Which parameter to set to the x axis in the QC plot (Default: nFeature_RNA).
#' @param y Which parameter to set to the y axis in the QC plot (Default: percent.mito).
#' @param log.scale Whether to plot with log.scale or not (Default: TRUE)
#' @param ncol Number of columns (Default = NA, determined from the data)
#' @param pt.size The point size, passed to ggplot (Default: 1)
#' @return A ggplot/cowplot grid object.
#' @keywords QC plot signatures highlight
#' @export
#' @examples
#' # After MySeuratObject has been log-normalized in order to contain an RNA>data assay>slot and been augmented with a percent.mito metadata column:
#' # SignatureList should be a list of character vectors, each containing a series of feature/gene names.
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject <- ScoreSignatures(MySeuratObject,SignatureList)
#' SignaturePlotQC(MySeuratObject,names(SignatureList))

SignaturePlotQC <- function(object, signatures, x= "nFeature_RNA", y="percent.mito", log.scale=TRUE, ncol=NA, pt.size=1){
  if(is.list(signatures)){
    sign.names <- names(signatures)
  }else if(is.vector(signatures, mode="character")){
    sign.names <- signatures
  }else{
    warning("The format of the 'signatures' argument (i.e. ",typeof(signatures),") is not supported. Provide a vector of metadata column names, or a named signature list.")
    return(ggplot())
  }
  sign.not.found <- setdiff(sign.names, colnames(object@meta.data))
  if(length(sign.not.found)>0){
    warning(paste("Signatures not found:",toString(sign.not.found)))
    sign.names <- intersect(sign.names, colnames(object@meta.data))
  }
  if(x %in% colnames(object@meta.data) & y %in% colnames(object@meta.data)){
    p <- list()
    for(n in sign.names){
      if(log.scale){
        p[[n]] <- ggplot(object@meta.data) + geom_point(aes_string(x=x, y=y, color=n), size=pt.size) + scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10') + lims(colour=c(0,NA))
      }else{
        p[[n]] <- ggplot(object@meta.data) + geom_point(aes_string(x=x, y=y, color=n), size=pt.size) + lims(colour=c(0,NA))
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
        warning(paste("Incorrect number of signatures found (",length(sign.name),")."))
      }
    }
    p <- cowplot::plot_grid(plotlist = p, ncol= ncol)
    return(p)
  }else{
    warning("Columns x or y (default: nFeature_RNA and percent.mito) not present in the object metadata.")
  }
}

