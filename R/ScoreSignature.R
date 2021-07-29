#' @export
ScoreSignature <- function(object, features, name = "Signature", assay=NA, slot.use="data"){
  print(paste0("Scoring the signature: ",name))
  if(is.na(assay)){
    if("imputed" %in% names(object@assays)){
      assay <- "imputed"
      print("Using the imputed assay.")
    }else if("RNA" %in% names(object@assays)){
      assay <- "RNA"
      print("Using the RNA assay.")
    }else{
      warning("No imputed or RNA assay found in the Seurat object. Provide an assay.")
    }
  }
  features.use <- intersect(features,rownames(slot(object[[assay]],slot.use)))
  missing.genes <- setdiff(unique(features), rownames(slot(object[[assay]],slot.use)))
  if(length(missing.genes)>0){
    print(paste("Missing features assumed null:",toString(missing.genes)))
  }
  scores <- colSums(rbind(slot(object[[assay]],slot.use)[features.use,],0))/length(features) # Concatenate a line of zeros because colSums doesn't deal with one-line matrices.
  object <- AddMetaData(object = object, metadata = scores, col.name = name)
}

ScoreSignatures <- function(object, signature.list, assay=NA, slot.use="data"){
  for(i in names(signature.list)){
    object <- ScoreSignature(object, signature.list[[i]], i, assay=assay, slot.use=slot.use)
  }
  return(object)
}
