#' WriteSignatures
#'
#' This function writes signature lists in the [.gmt format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29).
#' @param signatures Named list of chr vectors, in which the names of the list elements are the signature names, and the chr vectors are the corresponding gene collections.
#' @param file The address/name of the file from which signatures are read.
#' @keywords write signatures gmt .gmt gene matrix transpose
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

WriteSignatures <- function(signature.list,file,comment="cell type gene signatures"){
  for(i in 1:length(signature.list)){
    readr::write_tsv(as.data.frame(t(as.matrix(c(names(signature.list)[i],comment,signature.list[[i]])))),file,append = T,col_names = F)
  }
}
