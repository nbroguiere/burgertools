#' Load TCGA data
#'
#' Read a tsv file with TCGA most commonly mutated genes for a given cohort.
#'
#' Read a tab separated value (tsv) file with most commonly mutated genes for a given cancer cohort according to the cancer genome atlas (TCGA). The tsv file should have columns with Simple Somatic Mutation (SSM) Affected Cases in cohort, used to derive a percentage of patients affected, as well as gene (ensembl) IDs. This import function follows the conventions of files exported in tsv from:
#'
#' https://portal.gdc.cancer.gov/exploration?searchTableTab=genes
#'
#' @param file character(1), the tsv file path.
#' @return A data frame with at least gene_ensembl_id, cases_in_cohort and tcga_percent columns.
#' @keywords variants annotation TCGA common mutations cancer genome atlas
#' @examples
#' tcga_df <- ReadTCGA("Myfolder/MyTCGAcohortMostCommonMutatedGenes.tsv")
#' @export

ReadTCGA <- function(file){
  tcga <- readr::read_tsv(file,show_col_types = FALSE)
  colnames(tcga)[colnames(tcga)=="# SSM Affected Cases in Cohort"] <- "cases_in_cohort"
  colnames(tcga)[colnames(tcga)=="Gene ID"] <- "gene_ensembl_id"
  tcga$tcga_percent <- as.numeric(limma::strsplit2(x = tcga$cases_in_cohort,split = '\\(|%')[,2])
  tcga <- as.data.frame(tcga)
  rownames(tcga) <- tcga$gene_ensembl_id
  return(tcga)
}
