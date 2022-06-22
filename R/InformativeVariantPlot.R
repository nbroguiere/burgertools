#' Plot the informative variants with annotations
#'
#' Plots the excess entropy vs coverage of the variants, highlighting in color the informative variants. Optionally, can label all or only informative variants, or the variants having a minimal impact or present in a minimal percent of TCGA patients for the relevant cohort, or a combination thereof.
#'
#' @param object A seurat or genotype object on which FindInformativeVariants() has been computed.
#' @param do.label logical(1). Should variants be labeled (Default: TRUE).
#' @param label.column character(1). Which metadata column should be used as labels (Default: variants).
#' @param label.only.informative logical(1). Should the labels only be for informative variants (Default: TRUE).
#' @param min.impact character(1). Label only variants which have a predicted impact more than "LOW", "MODERATE", or "HIGH". NA to label all (Default:NA).
#' @param min.tcga.percent character(1). Label only variants which concern genes mutated in more than a certain percent of TCGA patients in the cohort used for genotype annotation. Requires a tcga_percent column to be present in the genotype object.
#' @param check.overlap logical(1). Should overlapping labels be omitted (Default: TRUE).
#' @param pt.size numeric(1). The point size (Default: 1).
#' @param text.size numeric(1). The text size (Default:3.5).
#' @param colors.use character(3). The colors for other variants, informative variants, and high impact / high tcga percent variants, in this order (Default: c("grey","#8888FF","red")).
#' @param log.scale logical(1). Should the axis be switched to log scale (Default: FALSE).
#' @param min.entropy numeric(1). Minimal excess_entropy value for variants to be represented by a dot in the graph. Values above 0, even small, greatly speed up rendering without noticeably changing the result, if many very low informative variants are present. More critical when using log scale, as low entropy variants are then visible (Default:1e-3).
#' @param assay character(1). If the object is a seurat object, which assay contains the consensus variant data (Default: "VAR").
#' @return Returns a ggplot
#' @keywords annotate variants vep genotype
#' @examples
#' InformativeVariantPlot(MyObject,label.column = "summary") # Other labels that can be interesting: symbol, summary, short_description.. see also colnames(MyGenotypes@metadata)
#' InformativeVariantPlot(MyObject,label.column = "summary",min.impact = "MODERATE") # Only label the variants that have a predicted impact more than moderate
#' InformativeVariantPlot(MyObject,label.column = "summary",min.impact = "MODERATE",min.tcga.percent = 1,label.only.informative=F,log.scale = T) # Add the restriction to only label above a minimal tcga percent. The results are few, so label even non-informative variants, and activate log scale to see better the low coverage variants.
#' @export
InformativeVariantPlot <- function(object, do.label=T, label.column="variants", label.only.informative=T, min.impact=NA, min.tcga.percent=NA, check.overlap=T, pt.size=1, text.size=3.5, colors.use=c("grey","#8888FF","red"),log.scale=F,min.entropy=1e-3,assay="VAR"){
  # If input is a Seurat object, convert to genotype object:
  if(class(object)=="Seurat"){
    object <- GenotypeObject(metadata = object[[assay]]@meta.features, variants = rownames(object[[assay]]), informative_variants = object[[assay]]@var.features)
  }

  # Apply min entropy
  object <- object[[object$excess_entropy>min.entropy,]]

  # Main dot plot
  other_variants <- setdiff(object@variants, object@informative_variants)
  if(label.column=="variants"){
    labels.tmp <- object@variants
  }else{
    labels.tmp <- object@metadata[,label.column,drop=T]
  }
  df <- data.frame(coverage=object$coverage, excess_entropy=object$excess_entropy, label=labels.tmp, row.names = object@variants)
  p <- ggplot(df,aes(x=coverage,y=excess_entropy,label=label))+
    geom_point(data=df[other_variants,],color=colors.use[1],size=pt.size)+
    geom_point(data=df[object@informative_variants,],color=colors.use[2],size=pt.size)

  # Handle the labels
  if(do.label){
    if(label.only.informative){
      variants_to_label <- object@informative_variants
    }else{
      variants_to_label <- object@variants
    }
    if(!is.na(min.impact)){
      variants_to_label <- intersect(variants_to_label,object@variants[object$impact>=min.impact])
    }
    if(!is.na(min.tcga.percent)){
      variants_to_label <- intersect(variants_to_label,object@variants[object$tcga_percent>min.tcga.percent])
    }
    p <- p + geom_text(data=df[variants_to_label,],size=text.size,check_overlap = check.overlap)
  }

  # Handle the "red" (or other third color) dots
  if(!is.na(min.impact) | !is.na(min.tcga.percent)){
    if(label.only.informative){
      variants_to_highlight <- object@informative_variants
    }else{
      variants_to_highlight <- object@variants
    }
    if(!is.na(min.impact)){
      variants_to_highlight <- intersect(variants_to_highlight,object@variants[object$impact>=min.impact])
    }
    if(!is.na(min.tcga.percent)){
      variants_to_highlight <- intersect(variants_to_highlight,object@variants[object$tcga_percent>min.tcga.percent])
    }
    p <- p + geom_point(data=df[variants_to_highlight,],color=colors.use[3],size=pt.size)
  }

  # Apply log scale
  if(log.scale){
    p <- p + scale_x_log10()+scale_y_log10()
  }
  return(p)
}
