% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Import10X.R
\name{Import10X}
\alias{Import10X}
\title{Import data in the 10X format with eventual metadata and dimensionality reductions into a Seurat object}
\usage{
Import10X(
  dir,
  project_name = "Project",
  assay = "RNA",
  gene.column = 1,
  cell.column = 1,
  unique.features = T,
  strip.suffix = T
)
}
\arguments{
\item{project_name}{character(1). A project name passed to CreateSeuratObject.Default: ""}

\item{assay}{character(1). An assay name passed to CreateSeuratObject. Default: "RNA"}

\item{gene.column}{integer(1). The column in the features.tsv that should be used as row names for the Seurat object. Passed to Read10X. Default: 1.}

\item{cell.column}{integer(1). The column in the barcodes.tsv that should be used as column names for the Seurat object. Passed to Read10X. Default: 1.}

\item{unique.features}{logical(1). Should the feature names be made unique. Passed to Read10X. Default: TRUE}

\item{strip.suffix}{logical(1). Should a constant suffix such as "-1" be stripped from the barcodes. Passed to Read10X. Default: TRUE}

\item{object}{Seurat object.}
}
\value{
A Seurat object.
}
\description{
Wrapper for Read10X that imports data in the 10X format (mtx raw count matrix, cells in columns, with features.tsv(.gz) and barcodes.tsv(.gz) files) directly into a Seurat object. Additionally looks for metadata.tsv(.gz) in the same folder and imports it as metadata of the Seurat object. If some columns of the form PC_i, UMAP_i, TSNE_i, CC_i, or HARMONY_i are found in the metadata, they are directly imported as dimensionality reductions (pca, umap, tsne, cca, and harmony respectively) in the Seurat object. The function combines well with Export10X in particular, to easily save and re-import 10X data that has been filtered using Seurat objects, and is helpful as well to share exported data in simple language-agnostic formats.
}
\examples{
# No project name, assay is RNA:
NewSeuratObject <- Import10X("MyDir")

# Re-export a Seurat object with two custom metadata columns, two dimensionality reductions, and without compression, then re-import it with a custom project and assay names:
Export10X(SeuratObject, "MyDir", c("nFeature_RNA", "mito.content"), c("pca","umap"), gzip=FALSE)
NewSeuratObject <- Import10X("MyDir","MyProject","mRNA3p")
}
\keyword{10X}
\keyword{Import}
\keyword{Reimport}
