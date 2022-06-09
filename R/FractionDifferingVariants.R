#' Fraction of differing variants
#'
#' Custom distance function tailored for variant data in a vartrix-like format (0=no data, 1=reference allele, 2=alternate allele, 3=heterozygous). The distance is defined as the fraction of variants which do not match, among the variants which are covered in both cells/genotypes being compared.
#'
#' @param x A (typically sparse) genotype matrix in vartrix conventions, with cells/samples as rows and variants as columns.
#' @return A distance object (as in stats package) containing the distances between the cells/samples.
#' @keywords fraction differing mutations variants distance SNPs genotypes vartrix
#' @export
#' @examples
#' # Typically used as a custom distance in combination with RunMDS in order to compute a dimensionality reduction based on variants:
#' DefaultAssay(MySeuratObject) <- "VAR" # variant calls, sparse matrix with values 1,2,3 for ref/alt/heterozygous in vartrix conventions. Assuming VariableFeatures have been set to most informative variants for this assay.
#' MySeuratObject <- RunMDS(MySeuratObject,FractionDifferingVariants)
#' DimPlot(MySeuratObject, reduction="mds")
#' MySeuratObject <- RunUMAP(MySeuratObject,reduction="mds",dims=1:30)
#' DimPlot(MySeuratObject, reduction="umap")
FractionDifferingVariants <- function(x){
  tmp <- cbind(x==1,x==2,x==3) # Could consider 2 and 3 (both "presence of variant") to be matching, instead of distinguishing heterozygous for variant / homozygous for variant.
  matching <- as.matrix(tmp %*%  t(tmp))
  x[x>0] <- 1
  overlap <- as.matrix(x %*% t(x))
  matching <- matching/overlap # Convert matching to a fraction instead of an absolute count, in place to save memory.
  matching[is.na(matching)] <- 0 # cells which had zero overlap are set to a shared mutation fraction of 0.
  return(stats::as.dist(1-matching))
}
