#' Gather features
#'
#' Gather into a dataframe several metadata columns and features scattered in various assays of a Seurat object. Accepts the format assay_feature to define a feature that is present in several assays non-ambiguously. Otherwise, will pick the features in priority from the default assay, or if not possible return an error.
#'
#' @param object Seurat object.
#' @param features character(n). Metadata columns and/or feature names, or assay feature pairs in the form assay_feature.
#' @return Returns a data frame with features or signatures / metadata columns as columns and cells as rows, in the same order as the cells/barcodes in the Seurat object.
#' @export
#' @examples
#' df <- GatherFeatures(MySeuratObject,c("nFeature_RNA","RNA_GAPDH"))

GatherFeatures <- function(object, features, slot="data"){
  # First pick up anything that can be picked up in the metadata:
  df <- object@meta.data[intersect(features,colnames(object@meta.data))]
  signatures_not_in_metadata <- setdiff(features,colnames(object@meta.data))

  # Then separate the features which have an assay defined (contain an underscore) from the rest:
  signatures_with_assay_defined <- grep("_",signatures_not_in_metadata,value=T)
  signatures_without_assay_defined <- setdiff(signatures_not_in_metadata,signatures_with_assay_defined)

  # Subdivide further the signatures without assay defined between those present in the DefaultAssay (picked in priority) and others:
  signatures_default_assay <- intersect(signatures_without_assay_defined,rownames(object[[DefaultAssay(object)]]))
  signatures_without_assay_defined <- setdiff(signatures_without_assay_defined,signatures_default_assay)

  # Gather the features in the default assay:
  if(length(signatures_default_assay)){
    df <- cbind(df,t(as.matrix(GetAssayData(object = object, assay = DefaultAssay(object), slot = slot)[signatures_default_assay,,drop=F])))
  }

  # Gather the features only given by name and not in the default assay:
  if(length(signatures_without_assay_defined)){
    for(j in 1:length(object@assays)){
      tmp <- intersect(signatures_without_assay_defined,rownames(GetAssayData(object = object, assay = Assays(object)[j], slot = slot)))
      if(length(tmp)){
        df <- cbind(df,t(as.matrix(GetAssayData(object = object, assay = Assays(object)[j], slot = slot)[tmp,,drop=F])))
      }
    }
  }

  # Gather the features given as assay_name:
  if(length(signatures_with_assay_defined)){
    for(i in 1:length(signatures_with_assay_defined)){
      tmp <- strsplit(signatures_with_assay_defined[i],"_")[[1]]
      defined.feat <- t(as.matrix(GetAssayData(object = object, assay = tmp[1], slot = slot)[tmp[2],,drop=F]))
      df <- cbind(df,defined.feat)
      colnames(df)[ncol(df)] <- signatures_with_assay_defined[i]
    }
  }

  # Check features are unique:
  if(length(colnames(df))!=length(unique(colnames(df)))){
    stop("Some signatures/features were not found were found in the metadata nor the default assay, and were present in more than 1 assay: \n", paste0(unique(colnames(df)[duplicated(colnames(df))]),collapse = " "),"\n",
         "Define them in a unique way declaring the assay from which they should be pulled with an underscore separator, as in: assay_feature. For example, RNA_GAPDH.")
  }

  return(df)
}
