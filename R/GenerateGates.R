#' Generate gates based on pair of thresholds (minimum main value and maximum of other values)
#'
#' Given a minimal value for the matching signature and maximum value for other signature, generate gates of the form: celltypeA = signatureA > min.main.value & signaturesNotA < max.other.value, celltypeB = signatureB > min.main.value & signaturesNotB < max.other.value, etc.
#'
#' The format of the gates returned is fit to be used directly in [ClassifyManual()], or can serve as a template to edit manually before using for [ClassifyManual()].
#'
#' @param signatures A list of signatures (in which case names of the list are used), or directly a vector of signature names.
#' @param min.main.value numeric(1). The minimal value for the signature matching the cell category. Default: 40.
#' @param max.other.value numeric(1). The maximal value for the signatures not matching the cell category. Default: 25.
#' @param population.names character(n). The name of the cell categories. Must be given in the same order as the name of the signatures with which they are associated. NA will use the signature names as population names.
#' @param display.gates logical(1). Should the gates generated (i.e. a string) be shown rather than returned silently. Default: TRUE.
#' @return Returns a string, character(1).
#' @keywords generate gates double threshold
#' @examples
#' MyGates <- GenerateGates(c("SignA","SignB","SignC"), min.main.value=1, max.other.value=0.1, population.names=c("CellA","CellB","CellC"))
#' @export

GenerateGates <- function(signatures, min.main.value=40, max.other.value=25, population.names=NA, display.gates=T){
  # Handle the situations in which signatures are the whole list of signatures, and we just use the names or is just a list instead of a vector but only names.
  if(is.list(signatures)){
    if(length(unlist(signatures))==length(signatures)){
      signatures <- unlist(signatures)
    }else{
      signatures <- names(signatures)
    }
  }else if(!is.vector(signatures, mode="character")){
    stop("The format of the 'signatures' argument (i.e. ",typeof(signatures),") is not supported. Provide a vector of metadata column names, or a list of names, or a named signature list.")
  }
  if(is.na(population.names[1])){
    population.names <- signatures
  }else{
    if(length(population.names)!=length(signatures)){
      stop("The number of signatures given does not match the number of population names.")
    }
  }
  gates <- list()
  for(i in 1:length(signatures)){
    gates[[i]] <- paste0(population.names[i]," = ",signatures[i],">",min.main.value, paste0("& ",signatures[-i],"<",max.other.value,collapse = " "))
  }
  gates_merged <- paste0(gates,collapse=",\n")
  if(display.gates){
    cat("Gates generated:\n")
    cat(gates_merged)
  }
  return(gates_merged)
}
