#' Scatter plots of signature scores on single cells
#'
#' Plots various signature scores against another signature score, in order to help define gates for a cell type. In the scatter plot, each point represents a cell in the Seurat object.
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

SignatureScatterPlot <- function(object,sign1,sign2,pt.size=0.8,n.rows=NA){
  if(length(sign1)>1 | !is.character(sign1)){
    stop("sign1 should be a character vector of length 1.")
  }
  sign2 <- unlist(sign2)
  n <- length(sign2)
  p <- list()
  for(i in 1:n){
    p[[i]] <- plotly::ggplotly(ggplot2::ggplot(cbind(object@meta.data,idents=Idents(object)))+ggplot2::geom_point(aes_string(sign1,sign2[i],color="idents"),size=pt.size)+ggplot2::ylab(label = sign2[i]))
    p[[i]] <- plotly::add_annotations(p[[i]], text = paste0(sign2[i]," vs ",sign1), x = 0.5, y = 1, yref = "paper", xref = "paper", xanchor = "middle", yanchor = "top", showarrow = FALSE, font = list(size = 15))

    if(i>1){
      p[[i]] <- plotly::style(p[[i]],showlegend=FALSE) # Make sure the legend is only shown once
    }
  }
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
