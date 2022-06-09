#' MAGIC imputation
#'
#' This function is a wrapper of the MAGIC imputation function to easily impute Seurat objects containing scRNA-seq data, with reasonable default parameters when used in the context of signature scoring for cell classification: only imputing the most variable genes and additional features of interest, and using only 40 PCs for efficiency, and using a narrow neighborhood to avoid over-smoothing.
#' @param object Seurat object. Must contain an RNA assay with a data slot that will be used for imputation.
#' @param features A subset of features/genes on which to run the imputation. Can be a list of signatures (Default: NULL/none, i.e. only impute all variable features if append.variable.features=TRUE as set by default).
#' @param append.variable.features Whether to append the variable features of the active assay in the Seurat object to the list of features to impute (Default: TRUE)
#' @param npca Number of principal components to use for imputation.
#' @param knn Number of nearest neighbors on which to compute bandiwth imputation.
#' @param t Diffusion parameter passed to MAGIC.
#' @param n.jobs Number of threads on which to run the imputation. Defaults to 6.
#' @param name Name of the assay storing the imputed data (Default: "imputed").
#' @param assay.use Name of the assay to impute (Default: "RNA").
#' @param slot.use Name of the slot to impute (Default: "data").
#' @return Seurat object with an additional slot containing MAGIC imputed data.
#' @keywords MAGIC Rmagic imputation
#' @export
#' @examples
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject[["imputed"]]

# Note: currently always using RNA and data assay, could extend to make it more flexible. Also give the choice of the name for the new slot where it ends up being stored.
Impute <- function(object, features=NULL, append.variable.features=TRUE, npca=40, knn=3, t=2, n.jobs=6, assay.use="RNA", slot.use="data", name="imputed"){
  if(!is.null(features)){
    features <- intersect(unique(unlist(features)), rownames(object))
    missing.genes <- setdiff(unique(unlist(features)), rownames(object))
    if(length(missing.genes)>0){
      warning(paste("The following genes were not found in the Seurat object and were omitted: \n",missing.genes))
    }
  }
  if(append.variable.features){
    if(length(VariableFeatures(object))){
      features <- unique(c(features,VariableFeatures(object)))
    }else{
      warning("Variable features required for imputation, but not present in the Seurat object. Imputation not computed. To run imputation on user-defined features only, set append.variable.features to FALSE.")
      return(object)
    }
  }
  if(length(features)>0){
    imputed <- Rmagic::magic(data = t(as.matrix(GetAssayData(object[features,], assay=assay.use, slot=slot.use))), npca=npca, knn=knn, t=t, n.jobs=n.jobs)
    imputed <- t(imputed$result)
    object[[name]] <- CreateAssayObject(imputed)
    return(object)
  }else{
    warning("None of the features required were found in the Seurat object. Imputation not computed.")
    return(object)
  }
}
