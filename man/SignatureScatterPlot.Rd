% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SignatureScatterPlot.R
\name{SignatureScatterPlot}
\alias{SignatureScatterPlot}
\title{Scatter plots of signature scores on single cells}
\usage{
SignatureScatterPlot(
  object,
  sign1,
  sign2,
  pt.size = 0.8,
  n.rows = NA,
  color.by = Idents(object)
)
}
\arguments{
\item{object}{Seurat object.}

\item{sign1}{character(1). The name of the signature to plot on the x axis.}

\item{sign2}{character(n). The name of the signature(s) to plot on the y axis.}

\item{pt.size}{numeric(1).The point size. Default: 0.8.}

\item{n.rows}{numeric(1).The number of rows in the plot grid. NA for automatic smart choice. Default: NA.}

\item{color.by}{character(n). The name of a metadata or feature column or a vector of idents or values to be used for point colors.}
}
\value{
A plotly grid of plots
}
\description{
Plots various signature scores (or features) against another signature score, in an interactive plot, in order to help define gates. In the scatter plot, each point represents a cell in the Seurat object. The points can be colored either by identity or according to a feature. Generalizes FeatureScatter() in that several scatter plots can be done at once, features and signatures can be pulled from several assays as well as from the metadata simultaneously, coloring can be done according to a feature instead of Idents(), and hovering can be used to get coordinates.
}
\examples{
SignatureNames <- names(MySignatureList)
SignatureScatterPlot(MySeuratObject,SignatureNames[1],SignatureNames[2])
SignatureScatterPlot(MySeuratObject,SignatureNames[1],SignatureNames[2:length(SignatureNames)])
SignatureScatterPlot(MySeuratObject,SignatureNames[1],SignatureNames[2:length(SignatureNames)],pt.size=1.2,n.rows=4)
}
\keyword{Plot}
\keyword{Scatter}
\keyword{Signature}
