#' MAGIC imputation
#'
#' This function is a wrapper of the MAGIC imputation function to easily impute Seurat objects containing scRNA-seq data.
#' @param object Seurat object. Must contain an RNA assay with a data slot that will be used for imputation.
#' @param features A subset of features on which to run the imputation. NA (default) uses VariableFeatures(object).
#' @param npca Number of principal components to use for imputation.
#' @param knn Number of nearest neighbors on which to perform data diffusion for imputation.
#' @param t Diffusion parameter passed to MAGIC.
#' @param n.jobs Number of threads on which to run the imputation. Defaults to 6.
#' @return Seurat object with an additional slot (accessed as object[["imputed"]]) containing MAGIC imputed data.
#' @keywords MAGIC Rmagic imputation
#' @export
#' @examples
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject[["imputed"]]

# Note: currently always using RNA and data assay, could extend to make it more flexible. Also give the choice of the name for the new slot where it ends up being stored.
Impute <- function(object, features = NA, npca=100, knn=3, t=2, n.jobs=6){
  if(is.na(features[1])){
    if(length(VariableFeatures(object))>0){
      features <- VariableFeatures(object)
    }
    else{
      object <- FindVariableFeatures(object)
      features <- VariableFeatures(object)
    }
  }else{
    features <- intersect(unique(features), rownames(object))
    missing.genes <- setdiff(unique(features), rownames(object))
    if(length(missing.genes)>0){
      warning(paste("The following genes were not found in the Seurat object and were omitted: \n",missing.genes))
    }
  }
  if(length(features)>0){
    imputed <- magic::magic(data = t(as.matrix(object@assays$RNA@data[features,])),npca=npca,knn=knn,t=t,n.jobs=n.jobs)
    imputed <- t(imputed$result)
    object[["imputed"]] <- CreateAssayObject(imputed)
    return(object)
  }else{
    warning("None of the features required were found in the Seurat object. Imputation not computed.")
    return(object)
  }
}
