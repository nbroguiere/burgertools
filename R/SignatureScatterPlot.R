#' Scatter plots of signature scores on single cells
#'
#' Plots various signature scores (or features) against another signature score, in an interactive plot, in order to help define gates. In the scatter plot, each point represents a cell in the Seurat object. The points can be colored either by identity or according to a feature. Generalizes FeatureScatter() in that several scatter plots can be done at once, features and signatures can be pulled from several assays as well as from the metadata simultaneously, coloring can be done according to a feature instead of Idents(), and hovering can be used to get coordinates.
#'
#' @param object Seurat object.
#' @param sign1 character(1). The name of the signature to plot on the x axis.
#' @param sign2 character(n). The name of the signature(s) to plot on the y axis.
#' @param pt.size numeric(1).The point size. Default: 0.8.
#' @param n.rows numeric(1).The number of rows in the plot grid. NA for automatic smart choice. Default: NA.
#'
#' @return A plotly grid of plots
#' @keywords Signature Scatter Plot
#' @export
#' @examples
#' SignatureNames <- names(MySignatureList)
#' SignatureScatterPlot(MySeuratObject,SignatureNames[1],SignatureNames[2])
#' SignatureScatterPlot(MySeuratObject,SignatureNames[1],SignatureNames[2:length(SignatureNames)])
#' SignatureScatterPlot(MySeuratObject,SignatureNames[1],SignatureNames[2:length(SignatureNames)],pt.size=1.2,n.rows=4)

SignatureScatterPlot <- function(object,sign1,sign2,pt.size=0.8,n.rows=NA,color.by=Idents(object)){
  if(length(sign1)>1 | !is.character(sign1[1])){
    stop("sign1 should be a character vector of length 1.")
  }
  sign2 <- unlist(sign2)
  n <- length(sign2)

  # If the colors are only a column name rather than a vector of cell idents, pick up the column in the metadata:
  if(length(color.by)==1){
    color.by <- object@meta.data[,color.by,drop=T]
  }

  # Gather the signatures and features to plot:
  df <- GatherFeatures(object,c(sign1,sign2))
  df <- cbind(df,color.by=color.by)

  # Plot
  p <- list()
  for(i in 1:n){
    p[[i]] <- plotly::ggplotly(ggplot2::ggplot(df)+ggplot2::geom_point(aes_string(paste0("`",sign1,"`"),paste0("`",sign2[i],"`"),color="color.by"),size=pt.size)+ggplot2::ylab(label = sign2[i]))
    p[[i]] <- plotly::add_annotations(p[[i]], text = paste0(sign2[i]," vs ",sign1), x = 0.5, y = 1, yref = "paper", xref = "paper", xanchor = "middle", yanchor = "top", showarrow = FALSE, font = list(size = 15))

    if(i>1){
      p[[i]] <- plotly::style(p[[i]],showlegend=FALSE) # Make sure the legend is only shown once
    }
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
  return(plotly::subplot(p, nrows = n.rows))
}
