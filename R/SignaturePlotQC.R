#' @export
SignaturePlotQC <- function(object, sign.names, x= "nFeature_RNA", y="percent.mito", log.scale=TRUE, ncol=NA){
  sign.not.found <- setdiff(sign.names, colnames(object@meta.data))
  if(length(sign.not.found)>0){
    warning(paste("Signatures not found:",toString(sign.not.found)))
    sign.names <- intersect(sign.names, colnames(object@meta.data))
  }
  if(x %in% colnames(object@meta.data) & y %in% colnames(object@meta.data)){
    p <- list()
    for(n in sign.names){
      if(log.scale){
        p[[n]] <- ggplot(object@meta.data) + geom_point(aes_string(x=x, y=y, color=n)) + scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10') + lims(colour=c(0,NA))
      }else{
        p[[n]] <- ggplot(object@meta.data) + geom_point(aes_string(x=x, y=y, color=n)) + lims(colour=c(0,NA))
      }
    }
    if(is.na(ncol)){
      if(length(sign.names)==1){
        ncol <- 1
      }else if(length(sign.names)>9){
        ncol <- 4
      }else if(length(sign.names)>4){
        ncol <- 3
      }else if(length(sign.names)>1){
        ncol <- 2
      }else{
        warning(paste("Incorrect number of signatures found (",length(sign.name),")."))
      }
    }
    p <- plot_grid(plotlist = p, ncol= ncol)
    return(p)
  }else{
    warning("Columns x or y (default: nFeature_RNA and percent.mito) not present in the object metadata.")
  }
}

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

IdentPlotQC <- function(object, ident=NA, x= "nFeature_RNA", y="percent.mito", log.scale=TRUE, ncol=NA){
  if(sum(is.na(ident))){
    object$tmp <- as.character(Idents(object))
    ident <- "tmp"
  }
  ident.not.found <- setdiff(ident, colnames(object@meta.data))
  if(length(ident.not.found)>0){
    warning(paste("Idents not found:",toString(ident.not.found)))
    ident <- intersect(ident, colnames(object@meta.data))
  }
  if(x %in% colnames(object@meta.data) & y %in% colnames(object@meta.data)){
    p <- list()
    for(n in ident){
      if(log.scale){
        p[[n]] <- ggplot(object@meta.data) + geom_point(aes_string(x="nFeature_RNA", y="percent.mito", color=n)) + scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')
      }else{
        p[[n]] <- ggplot(object@meta.data) + geom_point(aes_string(x="nFeature_RNA", y="percent.mito", color=n))
      }
    }
    if(is.na(ncol)){
      if(length(ident)==1){
        ncol <- 1
      }else if(length(ident)>9){
        ncol <- 4
      }else if(length(ident)>4){
        ncol <- 3
      }else if(length(ident)>1){
        ncol <- 2
      }else{
        warning(paste("Incorrect number of signatures found (",length(sign.name),")."))
      }
    }
    p <- plot_grid(plotlist = p, ncol= ncol)
    return(p)
  }else{
    warning("Columns x or y (default: nFeature_RNA and percent.mito) not present in the object metadata.")
  }
}
