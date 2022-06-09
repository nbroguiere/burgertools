#' Read a VEP tsv file
#'
#' This function reads a variant effect predictor (VEP) output file in tsv format, and returns the contents as a tibble.
#'
#' This function is a readr::read_tsv wrapper that reads a variant effect predictor (VEP) output file in tsv format, and returns the contents as a tibble. By default, removes uninformative columns (i.e. which have only one value), and setups the IMPACT as an ordered factor.
#'
#' @param file character(1) The path to the VEP file, can be compressed (.vcf or .vcf.gz).
#' @param remove_empty_col logical(1) Should uninformative columns with only a unique value/level be removed. Default: TRUE.
#' @return A tibble
#' @keywords genotypes variants vcf
#' @examples
#' genotypes <- ReadVEP("MyFolder/MyVEPresults.tsv.gz")
#' @export

ReadVEP <- function(file,remove_empty_col=TRUE){
  vep <- readr::read_tsv(file = file, comment = "##",show_col_types = FALSE)
  colnames(vep)[1] <- stringr::str_replace(colnames(vep)[1],pattern = "#",replacement = "")
  if(remove_empty_col){
    num_unique_el <- unlist(lapply(vep,function(x){return(length(unique(x)))}))
    vep[,num_unique_el<2] <- NULL
  }
  if("IMPACT" %in% colnames(vep)){
    vep$IMPACT <- factor(vep$IMPACT,levels=c("-","MODIFIER","LOW","MODERATE","HIGH"),ordered = T)
  }
  return(vep)
}
