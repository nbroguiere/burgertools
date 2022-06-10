#' Plot the informative variants with annotations
#'
#' Plots the excess entropy vs coverage of the variants, highlighting in color the informative variants. Optionally, can label all or only informative variants, or the variants having a minimal impact or present in a minimal percent of TCGA patients for the relevant cohort, or a combination thereof.
#'
#' @param genotype A genotype object
#' @param do.label logical(1). Should variants be labeled (Default: TRUE).
#' @param label.column character(1). Which metadata column in the genotype object should be used as labels (Default: variants).
#' @param label.only.informative logical(1). Should the labels only be for informative variants (Default: TRUE).
#' @param min.impact character(1). Label only variants which have a predicted impact more than "LOW", "MODERATE", or "HIGH". NA to label all (Default:NA).
#' @param min.tcga.percent character(1). Label only variants which concern genes mutated in more than a certain percent of TCGA patients in the cohort used for genotype annotation. Requires a tcga_percent column to be present in the genotype object.
#' @param check.overlap logical(1). Should overlapping labels be omitted (Default: TRUE).
#' @param pt.size numeric(1). The point size (Default: 1).
#' @param text.size numeric(1). The text size (Default:3.5).
#' @param colors.use character(3). The colors for other variants, informative variants, and high impact / high tcga percent variants, in this order (Default: c("grey","blue","red")).
#' @param log.scale logical(1). Should the axis be switched to log scale (Default: FALSE).
#' @param min.entropy numeric(1). Minimal excess_entropy value for variants to be represented by a dot in the graph. Values above 0, even small, greatly speed up rendering without noticeably changing the result, if many very low informative variants are present. More critical when using log scale, as low entropy variants are then visible (Default:1e-3).
#' @return Returns a ggplot
#' @keywords annotate variants vep genotype
#' @examples
#' InformativeVariantPlot(MyGenotypes,label.column = "summary") # Other labels that can be interesting: symbol, summary, short_description.. see also colnames(MyGenotypes@metadata)
#' InformativeVariantPlot(MyGenotypes,label.column = "summary",min.impact = "MODERATE") # Only label the variants that have a predicted impact more than moderate
#' InformativeVariantPlot(MyGenotypes,label.column = "summary",min.impact = "MODERATE",min.tcga.percent = 1,label.only.informative=F,log.scale = T) # Add the restriction to only label above a minimal tcga percent. The results are few, so label even non-informative variants, and activate log scale to see better the low coverage variants.
#' @export
InformativeVariantPlot <- function(genotype, do.label=T, label.column="variants", label.only.informative=T, min.impact=NA, min.tcga.percent=NA, check.overlap=T, pt.size=1, text.size=3.5, colors.use=c("grey","blue","red"),log.scale=F,min.entropy=1e-3){
  genotype <- genotype[[genotype$excess_entropy>min.entropy,T]]
  other_variants <- setdiff(genotype@variants, genotype@informative_variants)
  if(label.column=="variants"){
    labels.tmp <- genotype@variants
  }else{
    labels.tmp <- genotype@metadata[,label.column,drop=T]
  }
  df <- data.frame(coverage=genotype$coverage, excess_entropy=genotype$excess_entropy, label=labels.tmp, row.names = genotype@variants)
  p <- ggplot(df,aes(x=coverage,y=excess_entropy,label=label))+
    geom_point(data=df[other_variants,],color=colors.use[1],size=pt.size)+
    geom_point(data=df[genotype@informative_variants,],color=colors.use[2],size=pt.size)
  # Handle the labels
  if(do.label){
    if(label.only.informative){
      variants_to_label <- genotype@informative_variants
    }else{
      variants_to_label <- genotype@variants
    }
    if(!is.na(min.impact)){
      variants_to_label <- intersect(variants_to_label,genotype@variants[genotype$impact>=min.impact])
    }
    if(!is.na(min.tcga.percent)){
      variants_to_label <- intersect(variants_to_label,genotype@variants[genotype$tcga_percent>min_tcga_percent])
    }
    p <- p + geom_text(data=df[variants_to_label,],size=text.size,check_overlap = check.overlap)
  }
  # Handle the "red" (or other third color) dots
  if(!is.na(min.impact) | !is.na(min.tcga.percent)){
    if(label.only.informative){
      variants_to_highlight <- genotype@informative_variants
    }else{
      variants_to_highlight <- genotype@variants
    }
    if(!is.na(min.impact)){
      variants_to_highlight <- intersect(variants_to_highlight,genotype@variants[genotype$impact>=min.impact])
    }
    if(!is.na(min.tcga.percent)){
      variants_to_highlight <- intersect(variants_to_highlight,genotype@variants[genotype$tcga_percent>min_tcga_percent])
    }
    p <- p + geom_point(data=df[variants_to_highlight,],color=colors.use[3],size=pt.size)
  }
  if(log.scale){
    p <- p + scale_x_log10()+scale_y_log10()
  }
  return(p)
}
