#' @export
Classify <- function(object, sign.names, expected.values=NA, metadata.name="celltype"){
  if(sum(is.na(expected.values))){
    expected.values <- vector(mode="numeric",length=length(sign.names))+1
  }
  sign.not.found <- setdiff(sign.names, colnames(object@meta.data))
  if(length(sign.not.found)>0){
    warning(paste("Signatures not found:",toString(sign.not.found)))
    sign.names <- intersect(sign.names, colnames(object@meta.data))
  }
  if(length(expected.values)!=length(sign.names)){
    warning("The length of expected values does not match the length of the signature names. Aborting.")
    return(object)
  }
  if(length(intersect(names(expected.values),sign.names))!=length(sign.names)){
    if(length(intersect(names(expected.values),sign.names))==0){
      print("No matching names found in the expected values vector, assuming values are given in matching order.")
      names(expected.values) <- sign.names
    }else{
      warning("The names in the expected values vector do not match the ones in the signature names. Aborting.")
      return(object)
    }
  }
  sign.scores <- data.frame(row.names = rownames(object@meta.data))
  for(n in sign.names){
    sign.scores[,n] <- object@meta.data[,n]/expected.values[n]
  }
  object[[metadata.name]] <- colnames(sign.scores)[max.col(m = sign.scores)]
  Idents(object) <- object[[metadata.name]]
  return(object)
}
