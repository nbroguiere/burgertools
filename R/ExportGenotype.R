#' Reexport vartrix matrices from a genotype object
#'
#' Exports the reference and alternate count matrices, as well as the barcodes (colnames), variant names (rownames) and metadata/reference genotypes.
#'
#' @param genotype A genotype object
#' @param dir character(1). The directory in which the data should be exported. Created if non-existent.
#' @param meta_columns character(n). Metadata columns that should be exported from the genotype Object. Can also be "all", or empty vector c() for none. Default: "all".
#' @param append_genotypes logical(1). Whether the bulk genotypes should also be exported (appended to metadata). Default: TRUE.
#' @param gzip logical(1). Should the exported files be compressed. If compressed, the filenames are appended with the extension ".gz". Default: TRUE
#' @param rows character(1). Name of the exported file that lists row names. Default: vartrix_variants.tsv
#' @param cols character(1). Name of the exported file that lists column names. Default: vartrix_barcodes.tsv
#' @param ref_counts character(1). Name of the exported file that contains the reference allele sparse raw count matrix. Default: vartrix_ref_matrix.mtx
#' @param alt_counts character(1). Name of the exported file that contains the alternate allele sparse raw count matrix. Default: vartrix_alt_matrix.mtx
#' @param meta character(1). Name of the exported file that contains the metadata. Default: metadata.tsv
#' @keywords Reexport export genotype variants
#' @export
#' @examples
#' # By default keeping all the metadata columns, and no dimensionality reductions.
#' Export10X(SeuratObject, "MyDir")
#' # Include only two custom metadata columns, two dimensionality reductions, and do not compress the matrices.
#' Export10X(SeuratObject, "MyDir", c("nFeature_RNA", "mito.content"), c("pca","umap"), gzip=FALSE)

ExportGenotype <- function(genotype, dir, meta_columns = "all", append_genotypes = T, gzip=T, rows = "vartrix_variants.tsv", cols = "vartrix_barcodes.tsv", ref_counts = "vartrix_ref_matrix.mtx", alt_counts = "vartrix_alt_matrix.mtx", meta = "vartrix_metadata.tsv"){
  dir.backup <- getwd()
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)

  write(colnames(genotype@vartrix$REF), file = cols)
  write(rownames(genotype@vartrix$REF), file = rows)
  if(length(meta_columns)){
    if(meta_columns[1]=="all"){
      meta_columns <- colnames(genotype@metadata)
    }
  }
  meta.tmp <- genotype@metadata[,meta_columns,drop=F]
  if(append_genotypes){
    colnames(genotype@matrix) <- paste0("GTmatrix_",colnames(genotype@matrix))
    meta.tmp <- cbind(meta.tmp,genotype@matrix)
  }
  if(ncol(meta.tmp)){
    readr::write_tsv(meta.tmp,meta)
  }
  Matrix::writeMM(obj = genotype@vartrix$REF, file = ref_counts)
  Matrix::writeMM(obj = genotype@vartrix$ALT, file = alt_counts)
  if(gzip){
    R.utils::gzip(cols,  overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    R.utils::gzip(rows,  overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    R.utils::gzip(ref_counts,overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    R.utils::gzip(alt_counts,overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    if(ncol(meta.tmp)){
      R.utils::gzip(meta,  overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
    }
  }
  setwd(dir.backup)
}
