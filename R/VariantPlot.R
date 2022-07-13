#' Plots variants
#'
#' Plots variants as dot color on top of a dimensionality reduction used as x-y coordinates. By default, reference allele is blue, alternate allele is red, and heterozygous is purple.
#'
#' @param object Seurat object.
#' @param variants character(n). The name of the variant(s) to plot (as point color)
#' @param reduction character(1). The name of the dimensionality reduction that should be used as xy axis. If NA, will use umap, or tsne, or pca, in this order. Default: NA.
#' @param dims numeric(2). The dimensions within the dimensionality reduction that should be used as xy axis. Default: c(1,2).
#' @param pt.size numeric(1).The point size. Default: 1.
#' @param n.rows numeric(1).The number of rows in the plot grid. NA for automatic smart choice. Default: NA.
#' @param assay numeric(1). The assay from which the variant calls (in vartrix conventions) should be pulled. Default: "VAR".
#' @param colors character(n). Four colors, that will be used respectively for: no call, ref/ref, alt/alt, ref/alt calls.
#' @param order logical(1). Should the cells be ordered, with no call at the bottom and heterozygous on top, before plotting. Default: T.
#'
#' @return Returns a cowplot grid of plots
#' @keywords Variant plot FeaturePlot variants umap
#' @export
#' @examples
#' VariantPlot(MySeuratObject, "KRAS-12-G/D-MODERATE-chr12-25245350-C-T")
#' VariantPlot(MySeuratObject, VariableFeatures(MySeuratObject[["VAR"]])[1:9], reduction="mds", dims=c(1,3))

VariantPlot <- function(object, variants, reduction=NA, dims=c(1,2), pt.size=1, n.rows=NA, assay="VAR", colors=c("lightgrey","blue","red","purple"), order=T){

  DefaultAssay(object) <- assay

  # Gather the variants in a data frame and convert to factors
  if(is.na(variants[1]) | !length(variants)){
    stop("Incorrect or empty list of variants. Aborting.")
  }
  df <- GatherFeatures(object,variants)
  df <- as.data.frame(apply(df,2,as.factor))
  n <- ncol(df)
  if(is.na(variants[1]) | !length(variants)){
    stop("None of the variants have been found in this Seurat object. Aborting.")
  }


  # Gather the dimensionality reduction
  if(is.na(reduction)){
    if("umap" %in% Reductions(object)){
      reduction <- "umap"
    }else if("tsne" %in% Reductions(object)){
      reduction <- "tsne"
    }else if("mds" %in% Reductions(object)){
      reduction <- "mds"
    }else if("pca" %in% Reductions(object)){
      reduction <- "pca"
    }else{
      stop("No umap, tsne or pca dimensionality reduction found in the Seurat object. Please choose a reduction.")
    }
  }
  dr <- SO[[reduction]]@cell.embeddings[,dims]

  # Plot
  p <- list()
  for(i in 1:n){
    df2 <- cbind(df[,i,drop=F],dr)
    if(order){
      df2 <- df2[order(as.numeric(df[,i])),]
    }else{
      no.data.lines <- rownames(df2)[df2[,1,drop=T] == 0]
      ye.data.lines <- rownames(df2)[df2[,1,drop=T]  > 0]
      df2 <- df2[c(no.data.lines,sample(ye.data.lines)),]
    }
    p[[i]] <- ggplot2::ggplot(df2,aes_string(x = colnames(dr)[1], y = colnames(dr)[2], color = paste0("`",colnames(df)[i],"`"))) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values=setNames(colors,c("0","1","2","3"))) +
      Seurat::NoLegend() +
      ggplot2::ggtitle(label = colnames(df)[i])
  }

  # Decide the grid dimensions
  if(is.na(n.rows)){
    if(n<=2){
      n.rows <- 1
    }else if(n<=6){
      n.rows <- 2
    }else if(n<=12){
      n.rows <- 3
    }else if(n<=16){
      n.rows <- 4
    }else{
      n.rows <- floor(sqrt(n))
    }
  }

  return(cowplot::plot_grid(plotlist = p, nrow = n.rows))
}
