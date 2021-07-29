#' @export
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
