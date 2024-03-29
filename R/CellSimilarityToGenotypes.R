#' Compute the similarity of single cells in a Seurat object to the reference genotypes in a genotype object
#'
#' Computes the fraction of matching variants between single cells in a Seurat object and the reference genotypes in a genotype object. Store this similarity measure in a new assay within the Seurat object.
#'
#' The similarity is defined as the fraction of variants which have an identical call (i.e. ref/ref, alt/alt, or ref/alt) in individual cells compared to the reference genotype. Variants which are not covered by any read at the single cell level are ignored and do not contribute to the similarity value. Variants which are covered by only one read are always considered a match with a heterozygous reference genotype.
#'
#' @param seurat A seurat object
#' @param genotype A genotype object
#' @param assay.name character(1). The name under which the similarities are stored in the seurat object (Default: "similarity").
#' @param prefix character(1). A prefix appended to the genotype names to generate the feature names in the similarity assay (Default: "Similarity-").
#' @param features.use character(n). The features (i.e. variants) to be used for the similarity calculations. Can be a vector of variants by name, or NA to use informative variants if available, or "all" (Default: "all").
#' @param assays character(3). The name of the assays within the Seurat object that contain the variant calls (sparse matrix with values 1,2,3 in vartrix conventions), and counts of reference and alt alleles (Default: c("VAR","REF","ALT")).
#' @param layers The layers to use within the assays above (Default: c("data","data","data")).
#' @return Returns the Seurat object
#' @keywords compute cell similarity genotype fraction matching variants
#' @examples
#' MySeuratObject <- CellSimilarityToGenotypes(MySeuratObject,MyGenotypes)
#' DefaultAssay(MySeuratObject) <- "similarity"
#' FeaturePlot(MySeuratObject,"Similarity-MyGenotypeName1")
#' @export
CellSimilarityToGenotypes <- function(seurat, genotype, assay.name="similarity", prefix="Similarity-", features.use="all", assays=c("VAR","REF","ALT"), layers=c("data","data","data")){
  # Choose the features on which the similarity is computed (default: informative variants. If not present, all variants)
  if(is.na(features.use[1])){
    if(length(VariableFeatures(seurat[[assays[1]]]))){
      cat("Using the variable features from the Seurat object.\n")
      features.use <- VariableFeatures(seurat[[assays[1]]])
    }else if(length(genotype@informative_variants)){
      cat("Using the variable features from the Genotype object.\n")
      features.use <- genotype@informative_variants
    }else{
      warning("No informative variants / variable features found. Using all variants.\n")
      features.use <- genotype@variants
    }
  }
  if(features.use[1]=="all" & length(features.use)==1){
    cat("Using all variants present in the genotype object.\n")
    features.use <- genotype@variants
  }
  tmp <- setdiff(features.use,rownames(seurat[[assays[1]]]))
  if(length(tmp)){
    cat("Some variants present in the genotype object are missing from the Seurat object, excluding them.\n")
    cat("Number of missing variants:",length(tmp),"\n")
    features.use <- intersect(features.use,rownames(seurat[[assays[1]]]))
    cat("Number of variants remaining:",length(features.use),"\n")
  }
  tmp0 <- Matrix::t(GetAssayData(seurat, assay = assays[1], layer = layers[1])[features.use,])
  cvg0 <- Matrix::t(GetAssayData(seurat, assay = assays[2], layer = layers[2])[features.use,] + GetAssayData(seurat, assay = assays[3], layer = layers[3])[features.use,])
  tmp1 <- Matrix::t(genotype[features.use,])

  # # Previous simple version - perfect match only:
  # tmp2 <- cbind(tmp0==1,tmp0==2,tmp0==3)
  # tmp3 <- cbind(tmp1==1,tmp1==2,tmp1==3)

  # New improved version - when the count is 1, consider both the ref and the alt to be a match to a heterozygous reference genotype:
  tmp2 <- Matrix::drop0(cbind(tmp0==1,tmp0==2,tmp0==3 | cvg0==1))
  tmp3 <- Matrix::drop0(cbind(tmp1==1,tmp1==2,tmp1==3))

  matching <- as.matrix(tmp3 %*%  Matrix::t(tmp2))
  overlap <- as.matrix((tmp1>0) %*% Matrix::t(tmp0>0))
  matching <- Matrix::t(matching/overlap) # Convert matching to a fraction instead of an absolute count.
  matching[is.na(matching)] <- 0 # cells which had zero overlap are set to a shared mutation fraction of 0.
  colnames(matching) <- paste0(prefix, stringr::str_replace(colnames(matching),"_","-"))
  rownames(matching) <- colnames(seurat)

  cat("Storing similarities of single cells to reference genotypes in the assay '",assay.name,"'. \nExample of feature name: ",colnames(matching)[1], "\n", sep = "")
  seurat[[assay.name]] <- CreateAssayObject(counts = Matrix::t(matching))

  return(seurat)
}
