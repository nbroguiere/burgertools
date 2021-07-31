#' QC by cell type
#'
#' This function does cell-type aware quality control for scRNA-seq data in Seurat objects. Cells listed in 'celltypes' are filtered, keeping only cells for which 'parameter' is between min.val and max.val.
#'
#' @param object Seurat object. Must contain nFeature_RNA and percent.mito metadata columns, or the x and y plot parameters must be changed accordingly to plot other QC info otherwise.
#' @param celltypes Identities that should be filtered. Must be present in Idents(object).
#' @param parameter Metadata column on which to do the filtering.
#' @param min.val Minimal accepted values for each identity (Default: 0).
#' @param max.val Maximal accepted values for each identity (Default: +Inf).
#' @return A Seurat object
#' @keywords Filtering filter cell type aware celltype QC
#' @export
#' @examples
#' MySeuratObject <- RunQC(MySeuratObject)

RunQC <- function(object, celltypes, parameter, min.val=0, max.val=Inf){
  # Check the type of the inputs.
  # Check the inputs have the same length.
  # Check parameter is in the metadata, and the celltypes exist.
  for(i in 1:length(celltypes)){
    cells.to.exclude <- colnames(object)[(object[[parameter]]>max.val[i] | object[[parameter]]<min.val[i]) & Idents(object)==celltypes[i]]
    object <- object[,setdiff(colnames(object),cells.to.exclude)]
  }
  return(object)
}
