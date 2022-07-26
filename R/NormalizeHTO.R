#' Normalize hashtag oligo data
#'
#' Normalizes hashtag oligo (HTO) data. By default, first scales each HTO to its standard deviation to account for difference of intensity for the various labels, then scales each cell so that HTO are reported as a percentage. Alternatively, can scale the HTOs to a quantile of their values, or to their max after imputation.
#'
#' @param object Seurat object.
#' @param scale.HTO character(1) or logical(1). "stdev"/TRUE/T to scale each HTO to standard deviation across cells, No/no/FALSE/F for no scaling, "quantile" to scale to a quantile of values, and "imputed_max" to scale to the max after smoothing by imputation to avoid outliers. Default: "stdev".
#' @param normalize.cells numeric(1). The value to which each cell should be normalized. 0 for no cell-normalization. Default: 100.
#' @param assay character(1). The name of the assay which should be normalized. Default: "HTO".
#' @param q numeric(1). Only used if scale.HTO="quantile". The quantile used for normalization.
#' @param imputed.assay character(1). Only used if scale.HTO="imputed_max". The name of the assay that contains imputed HTO data.
#' @return Returns a Seurat object with normalized HTO data in the data slot.
#' @keywords HTO normalization
#' @export
#' @examples
#' MySeuratObject <- NormalizeHTO(MySeuratObject)

NormalizeHTO <- function(object, scale.HTO = "stdev", normalize.cells = 100, assay="HTO", q = 0.9, imputed.assay="imputed.HTO"){
  HTO_counts <- GetAssayData(object,"counts",assay)
  if(scale.HTO == "stdev" | scale.HTO==TRUE){
    HTO_counts <- sweep(HTO_counts, 1, matrixStats::rowSds(as.matrix(HTO_counts)), FUN = '/')   # Scale each row (HTO) to its stdev.
  }else if(scale.HTO == "quantile"){
    quantiles <- 1:nrow(HTO_counts)
    for(i in 1:nrow(HTO_counts)){
      quantiles[i] <- quantile(HTO_counts[i,], q)
    }
  }else if(scale.HTO == "imputed_max"){
    if(!imputed.assay %in% Assays(object)){
      stop(imputed.assay," is not present in the Seurat object.")
    }
    HTO_counts <- sweep(HTO_counts, 1, matrixStats::rowMaxs(as.matrix(object[[imputed.assay]]@data)), FUN = '/')   # Scale each row (HTO) to its max after imputation.
  }else if(scale.HTO == "no" | scale.HTO == "No" | scale.HTO==FALSE){
  }else{
    stop("Invalid value of scale.HTO. Valid options are stdev/TRUE/T, quantile, imputed_max, or No/no/FALSE/F.")
  }
  if(normalize.cells){
    HTO_counts <- sweep(HTO_counts, 2, matrixStats::colSums2(as.matrix(HTO_counts)), FUN = '/') # Scale each col (cbc) to its sum.
    HTO_counts <- HTO_counts*normalize.cells
  }
  object[["HTO"]]@data <- as(HTO_counts,"dgCMatrix")

  return(object)
}
