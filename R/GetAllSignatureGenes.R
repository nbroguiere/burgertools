#' Get all signature genes
#'
#' This function returns all the unique genes appearing within a list of signatures. Particularly useful when signatures contain "negative" genes, as it will remove the "-" sign. Also works on single vectors rather than lists, in which case it will remove the minus sign and duplicates. 
#' @param signature.list Character vector list. A list of signatures.
#' @return The unique elements present among the signatures (character vector). 
#' @keywords Unique genes signatures
#' @export
#' @examples
#' SL <- list(CD8TC=c("CD3","CD8","-CD4"), 
#'            Fibro=c("VIM","COL1A1","-EPCAM"))
#' UniqueGenes <- GetAllSignatureGenes(SL) # Returns c("CD3","CD4","CD8","VIM","COL1A1","EPCAM")
GetAllSignatureGenes <- function(signature.list){
  tmp1 <- grep("^-",unique(unlist(signature.list)),value=T,invert = T)
  tmp2 <- unlist(lapply(strsplit(grep("^-",unique(unlist(signature.list)),value=T),"^-"),`[`,2))
  tmp <- unique(c(tmp1,tmp2))
  return(tmp)
}
