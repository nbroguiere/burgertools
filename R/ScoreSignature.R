#' Signature scoring on Seurat objects
#'
#' This function scores a signature at the single cell level on Seurat objects (sum of the expression for each cell).
#' The signature must be a list of features (character vector).
#' Any feature can be set as a negative marker in the scoring by prefixing it with a minus sign, e.g. "-EPCAM". 
#' By default, the score is computed as the average expression of the features as normalized in the "data" slot (i.e. log-normalized typically). Can be overrun by setting the slot.use parameter. 
#' If imputed data is available (assay "imputed.RNA"), the imputed data will be used. Otherwise, data from the RNA assay will be used. This can be overrun by setting the assay.use parameter.
#' @param object Seurat object.
#' @param features The signature, typically a gene list (character vector). Negative features are prefixed by a minus sign as in "-EPCAM". 
#' @param name The name of the signature (character). The single cell scores are added to the Seurat object metadata in a column with this name, and can be accessed as object$name or displayed with FeaturePlot(object,name).
#' @param assay.use The assay within the Seurat object to use for signature scoring. If NA (default), it will use the "imputed" assay if it is present, or the "RNA" assay otherwise.
#' @param slot.use The slot to use within the assay to get data for signature scoring. By default, using the "data" assay.
#' @return Seurat object with an additional metadata slot containing the signature score.
#' @keywords Signature Scoring
#' @examples
#' #Scoring a single signature
#' MySeuratObject <- Impute(MySeuratObject)
#' MySeuratObject <- ScoreSignature(MySeuratObject,c("EPCAM","CDH1","ITGA6"),"EpithelialSignature")
#' head(MySeuratObject$EpithelialSignature)
#' FeaturePlot(MySeuratObject,"EpithelialSignature")
#' @export
#' 
ScoreSignature <- function(object, features, name = "Signature", assay.use=NA, slot.use="data"){
  print(paste0("Scoring the signature: ",name))
  if(is.na(assay.use)){
    if("imputed.RNA" %in% names(object@assays)){
      assay.use <- "imputed.RNA"
      print("Using the imputed.RNA assay.")
    }else if("RNA" %in% names(object@assays)){
      assay.use <- "RNA"
      print("Using the RNA assay.")
    }else{
      warning("No imputed.RNA or RNA assay found in the Seurat object. Explicitely provide an assay.")
    }
  }
  features.positive <- intersect(features,rownames(slot(object[[assay.use]],slot.use)))
  features.negative <- intersect(unlist(lapply(strsplit(grep("^-",features,value=T),"^-"),`[`,2)),rownames(slot(object[[assay.use]],slot.use)))
  
  missing.features <- setdiff(GetAllSignatureGenes(features), rownames(slot(object[[assay.use]],slot.use)))
  if(length(missing.features)>0){
    print(paste("Missing features assumed null:",toString(missing.features)))
  }
  scores.positive <- Matrix::colSums(rbind(slot(object[[assay.use]],slot.use)[features.positive,],0))/length(features) # Concatenate a line of zeros because colSums doesn't deal with one-line matrices.
  scores.negative <- Matrix::colSums(rbind(slot(object[[assay.use]],slot.use)[features.negative,],0))/length(features)
  
  object <- AddMetaData(object = object, metadata = scores.positive-scores.negative, col.name = name)
  return(object)
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
#' signature.list <- list(EpithelialSignature=c("EPCAM","CDH1","ITGA6"), FibroSignature=c("VIM","COL1A1","ITGA1","-EPCAM"))
#' MySeuratObject <- Impute(MySeuratObject,signature.list)
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
