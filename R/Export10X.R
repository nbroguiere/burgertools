#' Reexport the data from a Seurat object in 10X format
#'
#' Exports the counts matrix, features and barcodes from a Seurat object in a 10X-like format, with an additional metadata matrix in tsv format, that can be augmented with dimensionality reduction coordinates. The exported data is in standard mtx and tsv formats (sparse matrix matrix market format and tab separated values respectively), which facilitates sharing and reuse, and can be re-imported in a Seurat object with Import10X, or directly in matrices/data frames with e.g. Matrix::readMM and readr::read_tsv.
#'
#' @param object Seurat object.
#' @param dir character(1). The directory in which the data should be exported. Created if non-existent.
#' @param meta_columns character(n). Metadata columns that should be exported from the Seurat Object. Can also be "all", or empty vector c() for none. Default: "all".
#' @param append_reductions character(n). Name of the dimensionality reductions that should be included with the metadata. Default: "all".
#' @param gzip logical(1). Should the exported files be compressed. If compressed, the filenames are appended with the extension ".gz". Default: TRUE
#' @param rows character(1). Name of the exported file that lists row names. Default: features.tsv
#' @param cols character(1). Name of the exported file that lists column names. Default: barcodes.tsv
#' @param counts character(1). Name of the exported file that contains the sparse raw count matrix. Default: matrix.mtx
#' @param meta character(1). Name of the exported file that contains the metadata. Default: metadata.tsv
#' @param slot character(1). Name of the slot in the Seurat object from which the count matrix is exported. Default: counts
#' @param assay character(1). Name of the assay in the Seurat object from which the count matrix is exported. Default: RNA
#' @keywords Reexport export 10X
#' @export
#' @examples
#' # By default keeping all the metadata columns, and all dimensionality reductions.
#' Export10X(SeuratObject, "MyDir")
#' # Include only two custom metadata columns, two dimensionality reductions, and do not compress the matrices.
#' Export10X(SeuratObject, "MyDir", c("nFeature_RNA", "mito.content"), c("pca","umap"), gzip=FALSE)
#' # Export normalized CiteSeq data rather than raw RNA counts, and no metadata:
#' Export10X(SeuratObject, "MyDir", meta_columns = c(), slot="data", assay="CiteSeq")
Export10X <- function(object, dir, meta_columns = "all", append_reductions = "all", gzip=T, rows = "features.tsv", cols = "barcodes.tsv", counts = "matrix.mtx", meta = "metadata.tsv", slot="counts", assay="RNA"){
  dir.backup <- getwd()
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  data.tmp <- GetAssayData(object = object, slot = slot, assay = assay)
  write(colnames(data.tmp), file = cols)
  write(rownames(data.tmp), file = rows)
  if(length(meta_columns)){
    if(meta_columns[1]=="all"){
      meta_columns <- colnames(object@meta.data)
    }
  }
  meta.tmp <- object@meta.data[,meta_columns,drop=F]
  if(length(append_reductions)){
    if(append_reductions[1]=="all"){
      append_reductions <- Reductions(object)
    }
    for(i in 1:length(append_reductions)){
      meta.tmp <- cbind(meta.tmp,object[[append_reductions[i]]]@cell.embeddings)
    }
  }
  if(ncol(meta.tmp)){
    readr::write_tsv(meta.tmp,meta)
  }
  Matrix::writeMM(obj = data.tmp, file = counts)
  if(gzip){
    R.utils::gzip(cols,  overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    R.utils::gzip(rows,  overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    R.utils::gzip(counts,overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    if(ncol(meta.tmp)){
      R.utils::gzip(meta,  overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    }
  }
  setwd(dir.backup)
}
