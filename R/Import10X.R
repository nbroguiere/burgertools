#' Import data in the 10X format with eventual metadata and dimensionality reductions into a Seurat object
#'
#' Wrapper for Read10X that imports data in the 10X format (mtx raw count matrix, cells in columns, with features and barcodes tsv files, compressed or not) directly into a Seurat object. Additionally looks for a file named metadata.tsv or metadata.tsv.gz and imports it as metadata of the Seurat object. If some columns of the form PC_i, UMAP_i, TSNE_i, CC_i, or HARMONY_i are found in the metadata, they are directly imported as dimensionality reductions (pca, umap, tsne, cca, and harmony respectively) in the Seurat object. The function in particular combines well with Export10X, to easily save and re-import 10X data that has been filtered using Seurat objects, using simple, universal, reusable, sharable formats.
#'
#' @param object Seurat object.
#' @param project_name character(1). A project name passed to CreateSeuratObject.Default: ""
#' @param assay character(1). An assay name passed to CreateSeuratObject. Default: "RNA"
#' @param gene.column integer(1). The column in the features.tsv that should be used as row names for the Seurat object. Passed to Read10X. Default: 1.
#' @param cell.column integer(1). The column in the barcodes.tsv that should be used as column names for the Seurat object. Passed to Read10X. Default: 1.
#' @param unique.features logical(1). Should the feature names be made unique. Passed to Read10X. Default: TRUE
#' @param strip.suffix logical(1). Should the constant suffix such as -1 be stripped from the barcodes. Passed to Read10X. Default: TRUE
#' @return A Seurat object.
#' @keywords Reimport Import 10X
#' @export
#' @examples
#' # No project name, assay is RNA:
#' NewSeuratObject <- Import10X("MyDir")
#'
#' # Re-export a Seurat object with two custom metadata columns, two dimensionality reductions, and without compression, then re-import it with a custom project and assay names:
#' Export10X(SeuratObject, "MyDir", c("nFeature_RNA", "mito.content"), c("pca","umap"), gzip=FALSE)
#' NewSeuratObject <- Import10X("MyDir","MyProject","mRNA3p")

Import10X <- function(dir, project_name="", assay="RNA", gene.column=1, cell.column=1, unique.features=T, strip.suffix=T){
  dir_backup <- getwd()
  setwd(dir)

  counts <- Read10X(dir,gene.column=gene.column, cell.column=cell.column, unique.features = unique.features, strip.suffix = strip.suffix)

  if("metadata.tsv.gz" %in% list.files(dir) | "metadata.tsv" %in% list.files(dir)){
    meta_filename <- intersect(c("metadata.tsv.gz","metadata.tsv"),list.files(dir))[1]
    meta <- as.data.frame(readr::read_tsv(meta_filename,col_names = T))
    rownames(meta) <- colnames(counts)
    if("UMAP_1" %in% colnames(meta)){
      cat("\nFound UMAP in the metadata. Added to the Seurat object as a dimensionality reduction named umap.")
      ii <- grep("^UMAP_",colnames(meta))
      umap <- meta[,ii]
      meta <- meta[,-ii]
    }
    if("PC_1" %in% colnames(meta)){
      cat("\nFound PCA in the metadata. Added to the Seurat object as a dimensionality reduction named pca.")
      ii <- grep("^PC_",colnames(meta))
      pca <- meta[,ii]
      meta <- meta[,-ii]
    }
    if("CC_1" %in% colnames(meta)){
      cat("\nFound CC in the metadata. Added to the Seurat object as a dimensionality reduction named cca.")
      ii <- grep("^CC_",colnames(meta))
      pca <- meta[,ii]
      meta <- meta[,-ii]
    }
    if("HARMONY_1" %in% colnames(meta)){
      cat("\nFound Harmony in the metadata. Added to the Seurat object as a dimensionality reduction named harmony.")
      ii <- grep("^HARMONY_",colnames(meta))
      harmony <- meta[,ii]
      meta <- meta[,-ii]
    }
    if("TSNE_1" %in% colnames(meta)){
      cat("\nFound TSNE in the metadata. Added to the Seurat object as a dimensionality reduction named tsne")
      ii <- grep("^TSNE_",colnames(meta))
      tsne <- meta[,ii]
      meta <- meta[,-ii]
    }
    object <- CreateSeuratObject(counts = counts, project = project_name, assay = assay, min.cells = 0,min.features = 0,meta.data = meta)
    if(exists("umap"))
    {
      object[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), assay = assay)
    }
    if(exists("tsne"))
    {
      object[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne), assay = assay)
    }
    if(exists("cca"))
    {
      object[["cca"]] <- CreateDimReducObject(embeddings = as.matrix(cca), assay = assay)
    }
    if(exists("pca"))
    {
      object[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca), assay = assay)
    }
    if(exists("harmony"))
    {
      object[["harmony"]] <- CreateDimReducObject(embeddings = as.matrix(harmony), assay = assay)
    }
  }else{
    cat("Note that no metadata was found.")
    object <- CreateSeuratObject(counts = counts, project = project_name, assay = assay, min.cells = 0,min.features = 0)
  }

  setwd(dir_backup)
  return(object)
}
