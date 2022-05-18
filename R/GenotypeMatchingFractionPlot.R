GenotypeMatchingFractionPlot <- function(genotypes, key.title="Fraction of \nmatching mutations", limits=NULL, plot_method="plotly", colors = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"), ...){
  # Very important to handle the sparse matrices carefully to get it to run fast. Could get a further factor 2 by only computing half of the symmetrical matrix...
  tmp_matching <- as.matrix(t(genotypes@matrix==1) %*% (genotypes@matrix==1))+as.matrix(t(genotypes@matrix==2) %*% (genotypes@matrix==2))+as.matrix(t(genotypes@matrix==3) %*% (genotypes@matrix==3))
  tmp_overlap <- as.matrix(t(genotypes@matrix>0) %*% (genotypes@matrix>0))
  matching_fraction <- tmp_matching/tmp_overlap
  overlapping_fraction <- tmp_overlap/nrow(genotypes@matrix)
  return(heatmaply::heatmaply(matching_fraction, limits=limits, key.title=key.title, plot_method=plot_method, colors=colors, ...))
}
