#' @export
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
