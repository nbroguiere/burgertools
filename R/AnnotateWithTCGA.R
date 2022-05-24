#' Annotate with TCGA
#'
#' Add TCGA mutation prevalence within a cohort data to the metadata of a genotype object.
#'
#' See also ReadTCGA and https://portal.gdc.cancer.gov/exploration?searchTableTab=genes
#'
#' @param genotype A genotype object
#' @param tcga A data frame that must contain "gene_ensembl_id" and "tcga_percent" columns (character and numeric respectively).
#' @return A data frame with at least gene_ensembl_id, cases_in_cohort and tcga_percent columns.
#' @keywords vartrix variants matrix ref alt consensus freq frequency
#' @examples
#' MyGenotypes <- ReadVcf("MyGenotypes.vcf.gz")
#' MyTCGA <- ReadTCGA("MyTCGAmutations.tsv")
#' MyGenotypes <- AnnotateWithTCGA(MyGenotypes,MyTCGA)
#' @export

AnnotateWithTCGA <- function(genotype,tcga){
  genotype@metadata$tcga_percent <- 0
  rownames(tcga) <- tcga$gene_ensembl_id
  ii <- genotype@metadata$gene_ensembl_id %in% tcga$gene_ensembl_id
  genotype@metadata[ii,"tcga_percent"] <- tcga[genotype@metadata$gene_ensembl_id[ii],"tcga_percent",drop=T]
  return(genotype)
}
