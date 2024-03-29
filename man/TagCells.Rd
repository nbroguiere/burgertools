% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TagCells.R
\name{TagCells}
\alias{TagCells}
\title{Quickly tag cells in a Seurat Object}
\usage{
TagCells(
  object,
  plot.use = NULL,
  ident.name = "tag",
  ident.levels = c(TRUE, FALSE),
  ...
)
}
\arguments{
\item{object}{Seurat object. If using the default plot (DimPlot), a dimension reduction should have been computed, e.g. tsne or umap.}

\item{ident.name}{character(1) The name of the metadata column that will be filled the tag info.}

\item{ident.levels}{character(2). The names of the tagged vs other cells. Default: c(TRUE,FALSE)}

\item{...}{Other parameters passed to DimPlot when applicable.}
}
\value{
A Seurat object with an additional metadata column containing the cell type annotations and Idents() set to these annotations.
}
\description{
This function is a CellSelector wrapper that creates a new metadata column with two levels (TRUE and FALSE by default) and sets it as the new identity. Used as a convenience to track groups of cells across plots particularly during QC.
}
\examples{
SeuratObject <- TagCells(SeuratObject) # By default, tagging on a DimPlot.
SeuratObject[,SO$tag] # Seurat object restricted to the tagged cells
SeuratObject <- TagCells(SeuratObject,FeaturePlot(SeuratObject,"GAPDH")) # Can also pass a custom plot.
SeuratObject <- TagCells(SeuratObject,FeaturePlot(SeuratObject,"GAPDH"),ident.name="SelectedCells",ident.levels=c("Yes","no")) # Custom name for the new metadata column and levels.
}
\keyword{Cell}
\keyword{Classifier}
\keyword{celltype}
\keyword{classification}
\keyword{type}
