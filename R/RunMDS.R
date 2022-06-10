#' Run MDS
#'
#' Compute a multidimensional scaling (MDS) dimensionality reduction on a Seurat object, using the provided distance metric. By default, using 30 dimensions in the reduced space, based on the variable features of the current default assay and data slot, with results stored under the name "mds".
#'
#' Wrapper of the bigmds divide and conquer method, see [bigmds::divide_conquer_mds()] for details and references.
#'
#' For usage on variant data in vartrix-like sparse matrix conventions, use the custom distance function [FractionDifferingVariants()].
#'
#' For use on continuous rather than discrete data, typically consider [stats::dist()].
#'
#' @param object Seurat object.
#' @param distance Distance function, e.g. FractionDifferingVariants or stats::dist()
#' @param n.dims numeric(1). The number of dimensions in the reduction (Default:30).
#' @param features.use character(n). The features to use to compute the reduction. NA for VariableFeatures(object) (Default: NA).
#' @param assay character(1). The assay to use. NA for DefaultAssay (Default:NA).
#' @param slot character(1). The Seurat object assay slot to use (Default:"data").
#' @param reduction.name character(1). The name under which the dimensionality reduction will be stored under (Default: "mds").
#' @param key character(1). The key (prefix to column numbers) used to name the MDS coordinates (Default: "MDS_")
#' @param batch_size numeric(1). The number of cells processed for each batch of the divide and conquer MDS. The batches are then realigned and merged. Smaller batches provide faster runs, larger batches provide more accurate results (Default: 400).
#' @param n.cores numeric(1). The number of cores to use for parallel processing (Default: 1). Unix only.
#' @param seed numeric(1). Random seed set at the beginning of the run for reproducible results (Default:42).
#' @return Returns the Seurat object.
#' @keywords MDS multidimensional scaling dimensionality reduction dr Seurat wrapper
#' @export
#' @examples
#' DefaultAssay(MySeuratObject) <- "VAR" # variant calls, sparse matrix with values 1,2,3 for ref/alt/heterozygous in vartrix conventions. Assuming VariableFeatures have been set to most informative variants for this assay.
#' MySeuratObject <- RunMDS(MySeuratObject,FractionDifferingVariants)
#' DimPlot(MySeuratObject, reduction="mds")
#' MySeuratObject <- RunUMAP(MySeuratObject,reduction="mds",dims=1:30)
#' DimPlot(MySeuratObject, reduction="umap")
RunMDS <- function(object, distance, n.dims = 30, features.use = NA, assay=NA, slot="data", reduction.name="mds", key="MDS_", batch_size=400, n.cores=1, seed = 42){
  set.seed(seed)

  # Deal with default features and assays if not provided:
  if(is.na(features.use[1])){
    if(length(VariableFeatures(object))){
      features.use <- VariableFeatures(object)
    }else{
      warning("No features were provided, and no variable features were found in the current assay. Aborting MDS.")
      return(object)
    }
  }
  if(is.na(assay)){
    assay <- DefaultAssay(object)
  }

  # MDS with divide and conquer method:
  mds <- bigmds::divide_conquer_mds(x = t(GetAssayData(object, assay = assay, slot = slot)[informative_variants,]), l = batch_size, c_points = 2*n.dims, r = n.dims, n_cores = n.cores, dist_fn = distance)
  object[[reduction.name]] <- CreateDimReducObject(key = key, stdev = mds$eigen, assay = assay, embeddings = matrix(data = mds$points, nrow = nrow(mds$points), ncol = ncol(mds$points), dimnames = list(colnames(object),paste0(key,1:n.dims))))
  return(object)
}
