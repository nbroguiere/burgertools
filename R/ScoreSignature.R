#' Signature scoring on Seurat objects
#'
#' This function scores a signature at the single cell level on Seurat objects (sum of the expression for each cell).
#' The signature is a list of features (character vector).
#' By default, the score is computed as the average log-normalized expression of the features (i.e. using the "data" slot).
#' If imputed data is available (assay "imputed"), the imputed data will be used. Otherwise, data from the RNA assay will be used. This can be overrun by setting the assay parameter.
#' @param object Seurat object.
#' @param features The signature, typically a gene list (character vector).
#' @param name The name of the signature (character). The single cell scores are added to the Seurat object metadata in a column with this name, and can be accessed as object$name or displayed with FeaturePlot(object,name).
#' @param assay.use The assay within the Seurat object to use for signature scoring. If NA (default), it will use the "imputed" assay if it is present, or the "RNA" assay otherwise.
#' @param slot.use The slot to use within the assay to get data for signature scoring. By default, using the "data" assay.
#' @return Seurat object with an additional metadata slot containing the signature score.
#' @keywords Signature Scoring
#' @examples
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject <- ScoreSignature(MySeuratObject,c("EPCAM","CDH1","ITGA6"),"EpithelialSignature")
#' head(MySeuratObject$EpithelialSignature)
#' FeaturePlot(MySeuratObject,"EpithelialSignature")
#' @export

ScoreSignature <- function(object, features, name = "Signature", assay.use=NA, slot.use="data"){
  print(paste0("Scoring the signature: ",name))
  if(is.na(assay.use)){
    if("imputed" %in% names(object@assays)){
      assay.use <- "imputed"
      print("Using the imputed assay.")
    }else if("RNA" %in% names(object@assays)){
      assay.use <- "RNA"
      print("Using the RNA assay.")
    }else{
      warning("No imputed or RNA assay found in the Seurat object. Provide an assay.")
    }
  }
  features.use <- intersect(features,rownames(slot(object[[assay.use]],slot.use)))
  missing.genes <- setdiff(unique(features), rownames(slot(object[[assay.use]],slot.use)))
  if(length(missing.genes)>0){
    print(paste("Missing features assumed null:",toString(missing.genes)))
  }
  scores <- colSums(rbind(slot(object[[assay.use]],slot.use)[features.use,],0))/length(features) # Concatenate a line of zeros because colSums doesn't deal with one-line matrices.
  object <- AddMetaData(object = object, metadata = scores, col.name = name)
}

#' Multiple signature scoring on Seurat objects
#'
#' This function scores multiple signatures at the single cell level on a Seurat object (sum of the expression for each cell and signature).
#' By default, the score is computed as the average of the log-normalized expression of the features (i.e. using the "data" slot).
#' If imputed data is available (assay "imputed"), the imputed data will be used. Otherwise, data from the RNA assay will be used. This can be overrun by setting the assay paramter.
#' @param object Seurat object.
#' @param signature.list The named list of signatures. Each element of the list is a signature (character vector), the names of the elements in the list will be used as signature names.
#' @param assay.use The assay within the Seurat object to use for signature scoring. If NA (default), it will use the "imputed" assay if it is present, or the "RNA" assay otherwise.
#' @param slot.use The slot to use within the assay to get data for signature scoring. By default, using the "data" assay.
#' @return Seurat object with additional metadata slots containing signature scores.
#' @keywords Multiple Signature Scoring
#' @examples
#' MySeuratObject <- Impute(MySeuratObject)
#' signature.list <- list(c("EPCAM","CDH1","ITGA6"),c("VIM","COL1A1","ITGA1"))
#' names(signature.list) <- c("EpithelialSignature","FibroSignature")
#' MySeuratObject <- ScoreSignatures(MySeuratObject,signature.list)
#' head(MySeuratObject$EpithelialSignature)
#' head(MySeuratObject$FibroSignature)
#' FeaturePlot(MySeuratObject,c("EpithelialSignature","FibroSignature"))
#' @export

ScoreSignatures <- function(object, signature.list, assay.use=NA, slot.use="data"){
  for(i in names(signature.list)){
    object <- ScoreSignature(object, signature.list[[i]], i, assay.use=assay.use, slot.use=slot.use)
  }
  return(object)
}
