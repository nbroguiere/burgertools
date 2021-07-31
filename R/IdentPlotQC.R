#' QC plot with cell type highlight
#'
#' This function makes quality control (QC) plots from a Seurat object, by default percent.mito vs nFeature_RNA, while highlighting cell types in color.
#'
#' This is useful when making scRNA-seq from heterogeneous samples: for example immune cells are small and might be discarded as debris from much larger epithelial cancer cells if one is not aware of the differences in QC parameters occuring between various cell types.
#' @param object Seurat object. Must contain nFeature_RNA and percent.mito metadata columns, or the x and y plot parameters must be changed accordingly to plot other QC info otherwise.
#' @param ident The name of the categorical metadata column to use as cell types (Default: NA, using the current active annotations, i.e. Idents(object)). If ident is a list, an array of plots for the various set of annotations provided will be returned.
#' @param x Which parameter to set to the x axis in the QC plot (Default: nFeature_RNA).
#' @param y Which parameter to set to the y axis in the QC plot (Default: percent.mito).
#' @param log.scale Whether to plot with log.scale or not (Default: TRUE)
#' @param ncol Number of columns (Default = NA, determined from the data)
#' @param pt.size The point size, passed to ggplot (Default: 1)
#' @return A ggplot/cowplot grid object.
#' @keywords QC plot with cell types highlighted.
#' @export
#' @examples
#' # After MySeuratObject has been log-normalized in order to contain an RNA>data assay>slot and been augmented with a percent.mito metadata column:
#' # SignatureList should be a list of character vectors, each containing a series of feature/gene names. In this example, the signatures correspond to cell types, and the name of each signature is the name of the corresponding cell type (e.g. Tcell, Bcell, EpithelialCell, Fibroblasts, etc).
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject <- ScoreSignatures(MySeuratObject,SignatureList)
#' MySeuratObject <- Classify(MySeuratObject,names(SignatureList),"CellType") # Automatic cell type annotation based on cell type signatures.
#' IdentPlotQC(MySeuratObject,"celltype")

IdentPlotQC <- function(object, ident=NA, x= "nFeature_RNA", y="percent.mito", log.scale=TRUE, ncol=NA, pt.size=1){
  if(sum(is.na(ident))){
    object$tmp <- as.character(Idents(object))
    ident <- "tmp"
  }
  ident.not.found <- setdiff(ident, colnames(object@meta.data))
  if(length(ident.not.found)>0){
    warning(paste("Idents not found:",toString(ident.not.found)))
    ident <- intersect(ident, colnames(object@meta.data))
  }
  if(x %in% colnames(object@meta.data) & y %in% colnames(object@meta.data)){
    p <- list()
    for(n in ident){
      if(log.scale){
        p[[n]] <- ggplot(object@meta.data) + geom_point(aes_string(x="nFeature_RNA", y="percent.mito", color=n), size=pt.size) + scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')
      }else{
        p[[n]] <- ggplot(object@meta.data) + geom_point(aes_string(x="nFeature_RNA", y="percent.mito", color=n), size=pt.size)
      }
    }
    if(is.na(ncol)){
      if(length(ident)==1){
        ncol <- 1
      }else if(length(ident)>9){
        ncol <- 4
      }else if(length(ident)>4){
        ncol <- 3
      }else if(length(ident)>1){
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
