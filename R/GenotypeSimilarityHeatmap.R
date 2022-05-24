#' Genotype similarity heatmap
#'
#' Plot a heatmap of the similarity between the various genotypes in a genotype object.
#'
#' The similarity is defined as the fraction of matching variants, i.e. the number of identical variants divided by the number of variants which are simultaneously called in both genotypes being compared.
#'
#' @param genotypes genotype object, contains several genotypes to be compared on the heatmap.
#' @param key.title character(1), the legend/key. Default: "Fraction of matching mutations"
#' @param limits numeric(2), the limits of the color scale (typically c(0,1), which are the extremes values similarities could take), or NULL for min/max in current data.
#' @param plot_method "plotly" or "ggplot". Default "plotly"
#' @param colors a color palette, default viridis.
#' @param ... other graphic parameters passed on to heatmaply.
#' @return A heatmaply plot.
#' @keywords genotypes variants vcf
#' @examples
#' GenotypeSimilarityHeatmap(MyGenotypeObject)
#' @export

GenotypeSimilarityHeatmap <- function(genotypes, key.title="Fraction of \nmatching mutations", limits=NULL, plot_method="plotly", colors = viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"), ...){
  # Very important to handle the sparse matrices carefully to get it to run fast. Could get a further factor 2 by only computing half of the symmetrical matrix...
  tmp_matching <- as.matrix(Matrix::t(genotypes@matrix==1) %*% (genotypes@matrix==1))+as.matrix(Matrix::t(genotypes@matrix==2) %*% (genotypes@matrix==2))+as.matrix(Matrix::t(genotypes@matrix==3) %*% (genotypes@matrix==3))
  tmp_overlap <- as.matrix(Matrix::t(genotypes@matrix>0) %*% (genotypes@matrix>0))
  matching_fraction <- tmp_matching/tmp_overlap
  overlapping_fraction <- tmp_overlap/nrow(genotypes@matrix)
  return(heatmaply::heatmaply(matching_fraction, limits=limits, key.title=key.title, plot_method=plot_method, colors=colors, ...))
}
