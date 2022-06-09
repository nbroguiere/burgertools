#' Quickly tag cells in a Seurat Object
#'
#' This function is a CellSelector wrapper that creates a new metadata column with two levels (TRUE and FALSE by default) and sets it as the new identity. Used as a convenience to track groups of cells across plots particularly during QC.
#'
#' @param object Seurat object. If using the default plot (DimPlot), a dimension reduction should have been computed, e.g. tsne or umap.
#' @param ident.name character(1) The name of the metadata column that will be filled the tag info.
#' @param ident.levels character(2). The names of the tagged vs other cells. Default: c(TRUE,FALSE)
#' @param ... Other parameters passed to DimPlot when applicable.
#' @return A Seurat object with an additional metadata column containing the cell type annotations and Idents() set to these annotations.
#' @keywords Cell type classification celltype Classifier
#' @export
#' @examples
#' TagCells(SeuratObject) # By default, tagging on a DimPlot.
#' TagCells(SeuratObject,FeaturePlot(SeuratObject,"GAPDH")) # Can also pass a custom plot.
#' TagCells(SeuratObject,FeaturePlot(SeuratObject,"GAPDH"),ident.name="SelectedCells",ident.levels=c("Yes","no")) # Custom name for the new metadata column and levels.

TagCells <- function(object, plot="DimPlot", ident.name = "tag", ident.levels = c(TRUE,FALSE), ...){
  if(plot=="DimPlot"){
    cells.tmp <- Seurat::CellSelector(Seurat::DimPlot(object, ...))
  }else{
    cells.tmp <- Seurat::CellSelector(plot)
  }
  object@meta.data[,ident.name] <- ident.levels[2]
  object@meta.data[cells.tmp,ident.name] <- ident.levels[1]
  Idents(object) <- object@meta.data[,ident.name,drop=T]
  return(object)
}
