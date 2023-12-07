#' MAGIC imputation
#'
#' This function is a wrapper of the MAGIC imputation function to easily impute Seurat objects containing scRNA-seq data, with reasonable default parameters when used in the context of signature scoring for cell classification: only imputing the most variable genes and additional features of interest, and using only 40 PCs for efficiency, and using a narrow neighborhood to avoid over-smoothing.
#' @param object Seurat object. Must contain an RNA assay with a data layer (or a custom layer specified in argument layer.use) that will be used for imputation.
#' @param features A subset of features/genes on which to run the imputation. Can be a list of signatures in which the features in the list will be added, or can be "all" (Default: NULL/none, i.e. only impute all variable features if append.variable.features is TRUE, as set by default).
#' @param append.variable.features Whether to append the variable features of the active assay in the Seurat object to the list of features to impute (Default: TRUE)
#' @param npca Number of principal components to use for imputation (Default: 40).
#' @param knn Number of nearest neighbors on which to compute bandiwth imputation (Default: 3).
#' @param t Diffusion parameter passed to MAGIC (Default:2).
#' @param n.jobs Number of threads on which to run the imputation (Default: 6).
#' @param name Name of the assay storing the imputed data (Default: "imputed").
#' @param assay.use Name of the assay to impute (Default: "RNA").
#' @param layer.use Name of the layer to impute (Default: "data").
#' @param name Name of the assay in which the imputed data is stored (Default: paste0("imputed.",assay.use)).
#' @return Seurat object with an additional layer containing MAGIC imputed data.
#' @keywords MAGIC Rmagic imputation
#' @export
#' @examples
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject[["imputed.RNA"]]
#' SL <- list(CD8TC=c("CD3","CD8","-CD4"), 
#'            Fibro=c("VIM","COL1A1","-EPCAM"))
#' MySeuratObject <- Impute(MySeuratObject, features=SL, append.variable.features=TRUE, knn=4, t=3, npca=30, n.jobs=36, assay.use="RNA", layer.use="scale.data", name="custom.imputation")
#' MySeuratObject[["imputed.RNA"]]
Impute <- function(object, features=NULL, append.variable.features=TRUE, npca=40, knn=3, t=2, n.jobs=6, assay.use="RNA", layer.use="data", name=paste0("imputed.",assay.use)){

  if(!"Rmagic" %in% rownames(installed.packages())){
    stop('Rmagic is not installed, but needed for imputations. Consider running: devtools::install_github("cran/Rmagic"). The python version magic-impute is also needed, see install instructions at https://github.com/cran/Rmagic.')
  }

  if(!Rmagic::pymagic_is_available()){
    stop("Rmagic is installed, but pymagic is not available. Make sure Rmagic, reticulate, and pymagic are functional. Can be tested with Rmagic::pymagic_is_available().")
  }

  backup_default_assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay.use

  if(assay.use=="RNA" & DefaultAssay(object)!="RNA"){
    cat("Beware that the imputed assay is RNA, even though the default assay is currently ",DefaultAssay(object),"\n")
  }

  if(!is.null(features)){
    if(features[1]=="all"){
      cat("Imputing all features.\n")
      features <- rownames(object)
    }

    requested.features <- GetAllSignatureGenes(features)
    missing.features <- setdiff(requested.features, rownames(object))
    features <- intersect(requested.features, rownames(object))
    if(length(missing.features)>0){
      cat(paste("The following features were not found in the Seurat object and were omitted: \n",paste(missing.features,collapse=" "),"\n\n"))
    }
  }

  if(append.variable.features){
    if(length(VariableFeatures(object))){
      features <- unique(c(features,VariableFeatures(object)))
    }else{
      stop("Variable features required for imputation, but not present in the Seurat object. Imputation not computed. To run imputation on user-defined features only, set append.variable.features to FALSE.")
    }
  }

  if(length(features)>0){
    imputed <- Rmagic::magic(data = t(as.matrix(GetAssayData(object, assay=assay.use, layer=layer.use)[features,])), npca=npca, knn=knn, t=t, n.jobs=n.jobs)
    imputed <- t(imputed$result)
    object[[name]] <- CreateAssayObject(imputed)
    DefaultAssay(object) <- backup_default_assay
    return(object)
  }else{
    stop("None of the features required were found in the Seurat object. Imputation not computed.")
  }
}
