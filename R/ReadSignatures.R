#' Read signatures in .gmt format
#'
#' This function reads signatures in the [.gmt format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) and returns them as a named list.
#' @param file character(1). The address/name of the file from which signatures are read.
#' @param discard.comments logical(1).  Whether the comment column of the gmt file should be discarded. Default: TRUE.
#' @return Returns a mamed list of chr vectors, in which the names of the list elements are the signature names, and the chr vectors are the corresponding gene collections. If discard comments=FALSE, will return a list of two in which the first element are the signatures, and the second element the comments.
#' @keywords read signatures gmt .gmt gene matrix transpose
#' @export
#' @examples
#' signature.list <- list()
#' signature.list[["TC"]] <- c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","GZMA","IFNG","NCAM1","KLRK1","JAK3")
#' signature.list[["BC"]] <- c("IGHM","CD79A","CD79B","MZB1","DERL3","IGLL5","MS4A1","IGHD","IGHG1","IGHG3")
#' signature.list[["MP"]] <- c("CD74","HLA-DRA","HLA-DRB1","CD68","ITGAX","CD14","FCGR1A","MRC1","LYZ")
#' signature.list[["EP"]] <- c("EPCAM","PHGR1","TFF3","CKB","AGR2","PLCG2","CDX2","PIGR","PHGR1","CD24","LGALS4","GPX2","CEACAM6","CLDN3","KRT8","KRT19","KRT18","TSPAN8","OLFM4","FABP1","REG1A")
#' signature.list[["CAF"]]<- c("LOXL2","TAGLN","COL6A1","ITGA1","MMP1","MMP3","ACTA2","COL1A1","CALD1","COL3A1","COL1A2","COL6A3","COL6A2","THY1","DCN","COL5A2","COL5A1","PDGFRB")
#' signature.list[["EC"]] <- c("CD34","PECAM1","LYVE1","PLVAP","RAMP2","VWF","ESM1","FLT1","PODXL","GNG11","ENG","SLC9A3R2","HSPG2")
#' WriteSignatures(signature.list,"MyFile.gmt")
#' signature.list2 <- ReadSignatures("MyFile.gmt")
#' library(zeallot)
#' c(signature.list2,comments) %<-% ReadSignatures("MyFile.gmt",discard.comments=F)

ReadSignatures <- function(file, discard.comments){
  rawsign <- readLines(con = file)
  rawsign <- strsplit(rawsign,split = "\t")
  signnames <- list()
  sign <- list()
  comments <- list()
  for(i in 1:length(rawsign)){
    signnames[[i]] <- rawsign[[i]][1]
    comments[[i]] <- rawsign[[i]][2]
    sign[[i]] <- rawsign[[i]][3:length(rawsign[[i]])]
  }

  # If all comments can be converted to numeric without NA, return them as numeric.
  defaultW <- getOption("warn")
  options(warn = -1)
  if(!sum(is.na(as.numeric(comments)))){
    comments <- as.numeric(comments)
  }
  options(warn = defaultW)

  names(sign) <- signnames
  names(comments) <- signnames

  if(discard.comments){
    return(sign)
  }else{
    return(list(sign,unlist(comments)))
  }
}
