#' QC plot with cell type highlight
#'
#' This function makes interactive quality control (QC) plots from a Seurat object, by default percent.mito vs nFeature_RNA, while color-coding various cell types.
#'
#' Plotly interactions: click on legend categories to toggle the display on/off, double click to isolate one cell type/category/ident.
#' Hover to get cell info (useful to choose filtering parameters), and select pan/zoom in the top-right menu to do close-ups and navigate.
#'
#' This is useful when making scRNA-seq from heterogeneous samples: for example immune cells are small and might be discarded as debris from much larger epithelial cancer cells if one is not aware of the differences in QC parameters occuring between various cell types.
#' @param object Seurat object. Must contain nFeature_RNA and percent.mito metadata columns, or the x and y plot parameters must be changed accordingly to plot other QC info otherwise.
#' @param ident The name of the categorical metadata column to use as cell types (Default: NA, using the current active annotations, i.e. Idents(object)). If ident is a list, an array of plots for the various set of annotations provided will be returned.
#' @param x Which parameter to set to the x axis in the QC plot (Default: nFeature_RNA).
#' @param y Which parameter to set to the y axis in the QC plot (Default: percent.mito).
#' @param log.scale Whether to plot with log.scale or not (Default: TRUE)
#' @param ncol Number of columns (Default = NA, determined from the data)
#' @param pt.size The point size, passed to ggplot (Default: 1)
#' @param interactive Enable interactive plot? (Default: TRUE).
#' @param colors.use For interactive plots, a colorbrewer2.org palette name (e.g. "YlOrRd" or "Blues"), or a vector of colors to interpolate in hexadecimal "#RRGGBB" format, or a color interpolation function like colorRamp(). For non-interactive plots, a vector of colors passed to ggplot2 scale_color_manual. If the vector is named, the values will be matched based on names.
#' @param ... Other arguments passed to plot_ly in the case of interactive plots.
#' @return If one ident is plotted and interactive is enabled, returns interactive plot (plotly). If several, returns a ggplot grid (cowplot).
#' @keywords QC plot with cell types highlighted.
#' @export
#' @examples
#' # After MySeuratObject has been log-normalized in order to contain an RNA>data assay>slot and been augmented with a percent.mito metadata column:
#' # SignatureList should be a list of character vectors, each containing a series of feature/gene names. In this example, the signatures correspond to cell types, and the name of each signature is the name of the corresponding cell type (e.g. Tcell, Bcell, EpithelialCell, Fibroblasts, etc).
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject <- ScoreSignatures(MySeuratObject,SignatureList)
#' MySeuratObject <- Classify(MySeuratObject,names(SignatureList),"CellType") # Automatic cell type annotation based on cell type signatures.
#' IdentPlotQC(MySeuratObject,"celltype")

IdentPlotQC <- function(object, ident=NA, x= "nFeature_RNA", y="percent.mito", log.scale=TRUE, ncol=NA, pt.size=1, interactive=TRUE, colors.use=NULL, ...){
  if(sum(is.na(ident))){
    object$tmp <- as.character(Idents(object))
    ident <- "tmp"
  }
  ident.not.found <- setdiff(ident, colnames(object@meta.data))
  if(length(ident.not.found)>0){
    warning(paste("Idents not found:",toString(ident.not.found)))
    ident <- intersect(ident, colnames(object@meta.data))
  }
  if(length(ident)==1){
    if(length(unique(object@meta.data[,ident]))==2 & is.null(colors.use)){
      colors.use <- c("grey","#4444FF")
    }
    if(length(unique(object@meta.data[,ident]))==1 & is.null(colors.use)){
      colors.use <- c("#4444FF")
    }
  }
  if(x %in% colnames(object@meta.data)){
    if(y %in% colnames(object@meta.data)){
      if(length(ident)==0){
        warning("No valid ident found. Abort.")
        return()
      }else if(length(ident)==1 & interactive){
        if(log.scale){
          p <- plotly::layout(plotly::plot_ly(data = object@meta.data,type = "scatter", x=as.formula(paste0("~",x)), y=as.formula(paste0("~",y)), color=as.formula(paste0("~",ident)), mode="markers", marker = list(size=3.5*pt.size), colors = colors.use, ...), xaxis=list(type="log"), yaxis=list(type="log"))
        }else{
          p <- plotly::plot_ly(data = object@meta.data,type = "scatter", x=as.formula(paste0("~",x)), y=as.formula(paste0("~",y)), color=as.formula(paste0("~",ident)), mode="markers", marker = list(size=3.5*pt.size), colors=colors.use, ...)
        }
      }else{
        p <- list()
        for(n in ident){
          if(log.scale){
            p[[n]] <- ggplot2::ggplot(object@meta.data) + geom_point(aes_string(x=x, y=y, color=n), size=pt.size) + scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')
          }else{
            p[[n]] <- ggplot2::ggplot(object@meta.data) + geom_point(aes_string(x=x, y=y, color=n), size=pt.size)
          }
          if(!is.null(colors.use)){
            p[[n]] <- p[[n]] + ggplot2::scale_color_manual(values = colors.use)
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
        if(length(p)==1){
          p <- p[[1]]
        }else{
          p <- cowplot::plot_grid(plotlist = p, ncol= ncol)
        }
      }
      return(p)
    }else{
      warning("The variable set on the y axis(",y,") is not present as a column in the Seurat object metadata.")
      return()
    }
  }else{
    warning("The variable set on the x axis(",x,") is not present as a column in the Seurat object metadata.")
    return()
  }
}
