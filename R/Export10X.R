#' Reexport the data from a Seurat object in 10X format
#'
#' Exports the counts matrix, features and barcodes from a Seurat object in a 10X-like format, with an additional metadata matrix in tsv format, that can be augmented dimensionality reduction coordinates. The exported data is in standard mtx and tsv formats, which facilitates sharing and reuse, and can be re-imported in a Seurat object with Import10X.
#'
#' @param object Seurat object.
#' @param dir character(1). The directory in which the data should be exported. Created if non-existent.
#' @param meta_columns character(n). Metadata columns that should be exported from the Seurat Object. NA for all, empty vector c() for none. Default: NA.
#' @param append_reductions character(n). Name of the dimensionality reductions that should be included with the metadata. Default: none.
#' @param gzip logical(1). Should the exported files be compressed. If compressed, the filenames are appended with the extension ".gz". Default: TRUE
#' @param rows character(1). Name of the exported file that lists row names. Default: features.tsv
#' @param cols character(1). Name of the exported file that lists column names. Default: barcodes.tsv
#' @param counts character(1). Name of the exported file that contains the sparse raw count matrix. Default: matrix.mtx
#' @param meta character(1). Name of the exported file that contains the metadata. Default: metadata.tsv
#' @keywords Reexport export 10X
#' @export
#' @examples
#' # By default keeping all the metadata columns, and no dimensionality reductions.
#' Export10X(SeuratObject, "MyDir")
#' # Include only two custom metadata columns, two dimensionality reductions, and do not compress the matrices.
#' Export10X(SeuratObject, "MyDir", c("nFeature_RNA", "mito.content"), c("pca","umap"), gzip=FALSE)

Export10X <- function(object, dir, meta_columns = NA, append_reductions = c(), gzip=T, rows = "features.tsv", cols = "barcodes.tsv", counts = "matrix.mtx", meta = "metadata.tsv"){
  dir.backup <- getwd()
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  data.tmp <- GetAssayData(object = object, slot = "counts", assay = "RNA")
  write(colnames(data.tmp), file = cols)
  write(rownames(data.tmp), file = rows)
  if(is.na(meta_columns[1])){
    meta_columns <- colnames(object@meta.data)
  }
  meta.tmp <- object@meta.data[,meta_columns,drop=F]
  if(length(append_reductions)){
    for(i in 1:length(append_reductions)){
      meta.tmp <- cbind(meta.tmp,object[[append_reductions[i]]]@cell.embeddings)
    }
  }
  readr::write_tsv(meta.tmp,meta)
  Matrix::writeMM(obj = data.tmp, file = counts)
  if(gzip){
    R.utils::gzip(cols,  overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    R.utils::gzip(rows,  overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    R.utils::gzip(counts,overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    R.utils::gzip(meta,  overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
  }
  setwd(dir.backup)
}
