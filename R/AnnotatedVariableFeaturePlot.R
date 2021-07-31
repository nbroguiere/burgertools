#' Annotated variable features plot
#'
#' This function displays a variable feature plot with feature annotations.
#'
#' @param object Seurat object. Must contain nFeature_RNA and percent.mito metadata columns, or the x and y plot parameters must be changed accordingly to plot other QC info otherwise.
#' @param which.labels Number of most variable features to label, or list of variable features to label. If NA, label all the variable features (Default: NA).
#' @param pt.size Size of the points on the plot (Default: 1)
#' @param text.size Size of the feature annotations (Default: 2.5)
#' @param log Plot the x-axis on log scale
#' @param xnudge Offset in x on label position (Default: 0)
#' @param ynudge Offset in y on label position (Default: 0)
#' @param colors Colors to specify non-variable/variable status (Default: c("Black","Cyan")).
#' @param selection.method Which method to pull (refer to ?Seurat::VariableFeaturePlot()).
#' @param assay Assay to pull variable features from
#' @return A ggplot object.
#' @keywords Variable features genes plot annotations annotated QC
#' @export
#' @examples
#' SeuratObject <- FindVariableFeatures(SeuratObject)
#' # Annotate all variable features:
#' AnnotatedVariableFeaturePlot(SeuratObject)
#' # Annotate a few user-defined genes:
#' AnnotatedVariableFeaturePlot(SeuratObject,c("OLFM4","LGR5","CD44","GAPDH"))
#  # Annotate only a subset of the most variable genes:
#' AnnotatedVariableFeaturePlot(SeuratObject,500)

AnnotatedVariableFeaturePlot <- function(object, which.labels=NA, pt.size=1, text.size=2.5, log=NULL, xnudge=0, ynudge=0, colors=c("Black","Cyan"),selection.method=NULL, assay=NULL){
  if(sum(is.na(which.labels))){
    points <- VariableFeatures(object,assay = assay)
  }else if(is.vector(which.labels,mode="numeric")){
    points <- VariableFeatures(object,assay = assay)[1:min(which.labels[1],length(VariableFeatures(object,assay = assay)))]
  }else if(is.vector(which.labels,mode="character")){
    points <- intersect(which.labels,rownames(object))
  }
  LabelPoints(plot = VariableFeaturePlot(object,cols = colors, pt.size = pt.size, log = log, selection.method = selection.method, assay=assay), points = points, repel = FALSE, xnudge = 0, ynudge = 0,size=text.size)
}

