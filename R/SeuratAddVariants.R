#' Add the variant information contained in a genotype object to a Seurat object
#'
#' The genotype object should contain a consensus call (vartrix conventions), reference allele count, alternate allele count, and a frequency matrix (typically offset by one to distinguish 0 frequency from no data in an efficient sparse matrix format), as well as metadata about the variants. By default, the single cell matrices will be stored in assays named VAR, REF, ALT and FREQ respectively, and the metadata will be added to the VAR assay feature metadata.
#'
#' @param seurat Seurat object.
#' @param genotype Genotype object.
#' @param consensus character(1). The name of the Seurat assay that will store the consensus matrix (values 0 = no data, 1 = ref/ref, 2=alt/alt, 3=ref/alt).
#' @param reference character(1). The name of the Seurat assay that will store the reference allele count matrix.
#' @param alternate character(1). The name of the Seurat assay that will store the alternate allele count matrix.
#' @param frequency character(1). The name of the Seurat assay that will store the frequency matrix (0 = no data, 1-2 = ref only to alt only).
#' @keywords variant matrices vartrix seurat
#' @export
#' @examples
#' # MySeuratObject <- SeuratAddVariants(MySeuratObject, MyGenotypeObject)

SeuratAddVariants <- function(seurat, genotype, consensus="VAR", reference="REF", alternate="ALT", frequency="FREQ"){
  seurat[[consensus]] <- CreateAssayObject(counts = genotype@vartrix[["VAR"]][ ,colnames(seurat)], min.cells = 0)
  seurat[[reference]] <- CreateAssayObject(counts = genotype@vartrix[["REF"]][ ,colnames(seurat)], min.cells = 0)
  seurat[[alternate]] <- CreateAssayObject(counts = genotype@vartrix[["ALT"]][ ,colnames(seurat)], min.cells = 0)
  seurat[[frequency]] <- CreateAssayObject(counts = genotype@vartrix[["FREQ"]][,colnames(seurat)], min.cells = 0)
  seurat[[consensus]]@meta.features <- genotype@metadata
  return(seurat)
}
