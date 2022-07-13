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
#' @param distance Distance function, e.g. FractionDifferingVariants for variants or stats::dist() for common distances applying to continuous data. NA will use FractionDifferingVariants if the assay used is "VAR", euclidean dist otherwise.
#' @param n.dims numeric(1). The number of dimensions in the reduction (Default:30).
#' @param features.use character(n). The features to use to compute the reduction. NA for most variable features of assay (Default: NA).
#' @param assay character(1). The assay to use. NA uses "VAR" if available, otherwise the DefaultAssay (Default:NA).
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
#' DefaultAssay(MySeuratObject) <- "VAR" # variant calls, sparse matrix with values 1,2,3 for ref/alt/heterozygous in vartrix conventions. Assuming VariableFeatures have been set to most informative variants for this assay. Will also use VAR by default if assay=NA, so no need to explicitely change the default assay.
#' MySeuratObject <- RunMDS(MySeuratObject)
#' DimPlot(MySeuratObject, reduction="mds")
#' MySeuratObject <- RunUMAP(MySeuratObject,reduction="mds",dims=1:30)
#' DimPlot(MySeuratObject, reduction="umap")
RunMDS <- function(object, distance=NA, n.dims = 30, features.use = NA, assay=NA, slot="data", reduction.name="mds", key="MDS_", batch_size=400, n.cores=1, seed = 42){
  set.seed(seed)

  # Resolve the assay to use:
  if(is.na(assay)){
    if("VAR" %in% Assays(object)){
      assay <- "VAR"
    }else{
      assay <- DefaultAssay(object)
    }
  }

  # Resolve the distance to use:
  if(is.na(distance)){
    if(assay=="VAR"){
      distance <- burgertools::FractionDifferingVariants
    }else{
      distance <- stats::dist
    }
  }

  # Resolve the features to use:
  if(is.na(features.use[1])){
    if(length(VariableFeatures(object,assay=assay))){
      features.use <- VariableFeatures(object,assay=assay)
    }else{
      warning("No features were provided, and no variable features were found in the current assay. Aborting MDS.")
      return(object)
    }
  }

  # There's a bug if the batch size is close to the matrix size or smaller. Give a warning and choose a value that works.
  if(batch_size>ncol(object)/2){
    warning("The batch size is bigger than half the total number of cells in the Seurat object. Proceeding with one batch only.")
    batch_size <- ncol(object)
  }

  # MDS with divide and conquer method:
  mds <- bigmds::divide_conquer_mds(x = Matrix::t(GetAssayData(object, assay = assay, slot = slot)[features.use,]), l = batch_size, c_points = 2*n.dims, r = n.dims, n_cores = n.cores, dist_fn = distance)
  object[[reduction.name]] <- CreateDimReducObject(key = key, stdev = mds$eigen, assay = assay, embeddings = matrix(data = mds$points, nrow = nrow(mds$points), ncol = ncol(mds$points), dimnames = list(colnames(object),paste0(key,1:n.dims))))
  return(object)
}
