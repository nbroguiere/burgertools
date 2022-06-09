#' MAGIC imputation of features of an assay based on the distances in another assay
#'
#' This function imputes features in an assay based on (the merger with) another assay. Typically using the high resolution in the RNA data to impute more noisy additional assays.
#'
#' @param object Seurat object.
#' @param impute.assay character(1). The name of the assay which should be imputed, or "metadata", or "default_assay" (Default: "default_assay").
#' @param name character(1). The name under which the resulting imputed assay is stored in the Seurat object (Default: "crossimputed")
#' @param impute.features character(n). Features on which to run the imputation in the impute.assay, or "all" for all features in the impute.assay, or "variable_features" for the variable features in the impute.assay (Default:"all").
#' @param reference.assay character(1). The name of the assay which is appended to the impute.assay to guide imputation distances.
#' @param reference.features character(n). Features from the reference.assay to use to guide the imputation of the impute.assay, or "all" for all features in the reference.assay, or "variable_features" for the variable features in the reference.assay (Default:"variable_features").
#' @param npca integer(1). Number of principal components to use for imputation.
#' @param knn integer(1). Number of nearest neighbors on which to compute bandiwth imputation.
#' @param t integer(1). Diffusion parameter passed to MAGIC.
#' @param n.jobs integer(1). Number of threads on which to run the imputation (Default: 6).
#' @param slot.use character(1). Name of the slot to impute (Default: "data").
#' @return Seurat object with an additional slot containing MAGIC cross-imputed data.
#' @keywords MAGIC Rmagic imputation crossimputation cross-imputation
#' @export
#' @examples
#' MySeuratObject <- CrossImpute(MySeuratObject,"other_assay_to_impute_based_on_RNA_distances")
#' MySeuratObject[["crossimputed"]] # Check result

# Note: currently always using RNA and data assay, could extend to make it more flexible. Also give the choice of the name for the new slot where it ends up being stored.
CrossImpute <- function(object, impute.assay="default_assay", name="crossimputed", impute.features="all", reference.assay="RNA", reference.features="variable_features", npca=40, knn=3, t=2, n.jobs=6, slot.use="data"){

  backup.default.assay <- DefaultAssay(object)
  if(impute.assay=="default_assay"){
    impute.assay <- DefaultAssay(object)
  }

  if(impute.assay=="metadata"){
    if(sum(impute.features=="all")){
      warning("Define explicit metadata columns to impute, all is not accepted for metadata imputation.")
      return(object)
    }else if(sum(impute.features=="variable_features")){
      warning("Define explicit metadata columns to impute, variable_features is not accepted for metadata imputation.")
      return(object)
    }else{
      missing.features <- setdiff(unique(unlist(impute.features)), colnames(object@metadata))
      impute.features <- intersect(unique(unlist(impute.features)), colnames(object@metadata))
      if(length(missing.features)>0){
        warning(paste("The following features were not found in the impute.assay of the Seurat object and were omitted: \n",paste(missing.features,collapse=" ")))
      }
    }
    df <- object@metadata[,impute.features,drop=F]
  }else{
    DefaultAssay(object) <- impute.assay
    if(sum(impute.features=="all")){
      impute.features <- rownames(object)
    }else if(sum(impute.features=="variable_features")){
      impute.features <- VariableFeatures(object)
    }else{
      missing.features <- setdiff(unique(unlist(impute.features)), rownames(object))
      impute.features <- intersect(unique(unlist(impute.features)), rownames(object))
      if(length(missing.features)>0){
        warning(paste("The following features were not found in the impute.assay of the Seurat object and were omitted: \n",paste(missing.features,collapse=" ")))
      }
    }
    if(!length(impute.features)){
      warning("No features to impute available, aborting.")
      return(object)
    }else{
      df <- t(GetAssayData(object = object, slot = slot.use)[impute.features,,drop=F])
    }
  }

  DefaultAssay(object) <- reference.assay
  if(sum(reference.features=="all")){
    reference.features <- rownames(object)
  }else if(sum(reference.features=="variable_features")){
    reference.features <- VariableFeatures(object)
  }else{
    missing.features <- setdiff(unique(unlist(reference.features)), rownames(object))
    reference.features <- intersect(unique(unlist(reference.features)), rownames(object))
    if(length(missing.features)>0){
      warning(paste("The following features were not found in the reference.assay of the Seurat object and were omitted: \n",paste(missing.features,collapse=" ")))
    }
    if(!length(reference.features)){
      warning("No reference features available, aborting.")
      return(object)
    }else{
      df <- cbind(df,t(GetAssayData(object = object, slot = slot.use)[reference.features,,drop=F]))
    }
  }

  imputed <- Rmagic::magic(data = as.matrix(df), npca=npca, knn=knn, t=t, n.jobs=n.jobs)
  imputed <- t(imputed$result)
  object[[name]] <- CreateAssayObject(imputed)

  DefaultAssay(object) <- backup.default.assay
  return(object)
}
