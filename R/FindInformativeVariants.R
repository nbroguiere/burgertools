#' Find informative variants
#'
#' Computes the fraction of detection of variants at the single cell level, and entropy associated with the variants, based on a genotype and a Seurat objects. Stores this information in the genotype object metadata, and compute the most informative variants. Returns the updated genotype object.
#'
#' FindInformativeVariants populates the following metadata columns and slots in the genotype object:
#'
#'  $excess_entropy (see below),
#'
#'  $coverage (number of cells with data for the variant),
#'
#'  $coverage_frac (fraction of the cell with data for the variant),
#'
#'  @variants_by_coverage (variants sorted from top to least coverage),
#'
#'  @variants_by_information (variants sorted from top to least excess entropy),
#'
#'  @informative_variants (most informative variants, to be used similarly to VariableFeatures for continuous data within Seurat Objects for further analysis).
#'
#'  The amount of information provided by a variant at the single cell level is computed as the Shannon entropy, sum(-p*log2(p)), minus the minimal entropy contributed by the mere coverage of the variant, i.e. the entropy value if the data had just the two levels nodata and data. This is analogous to the excess variance for continuous data, translated to discrete variant data. This quantity is referred to as excess_entropy in the genotype metadata.
#'
#' @param genotype A genotype object.
#' @param seurat A seurat object.
#' @param n.variants numeric(1). The number of most informative variants stored in the @informative_variants slot (Default: 10000).
#' @param assay character(1). The name of the assay which stores the single cell variant calls data within the Seurat object, in vartrix formatting (sparse matrix with values 1-2-3 for ref-alt-het).
#' @return Returns the updated genotype object
#' @keywords most informative variants genotype variable features excess entropy
#' @examples
#' MyGenotype <- FindInformativeVariants(MyGenotype, MySeuratObject)
#' @export

FindInformativeVariants <- function(genotype, seurat, n.variants=10000, assay="VAR"){

  # Compute coverage
  vars_matrix <- seurat[[assay]]@counts
  coverage <- setNames(sparseMatrixStats::rowSums2(vars_matrix>0),rownames(vars_matrix))

  # # Apply a min.fraction (if including it as an additional function parameter, removed at the moment as the code runs very fast already anyway):
  # filter <- coverage >= ncol(seurat)*min.fraction
  # vars_matrix <- vars_matrix[filter,]
  # coverage <- coverage[filter]

  # Compute entropy associated with mutations as a way to discover the most informative ones:
  ni <- list()
  pi <- list()
  log2pi <- list()
  for(i in 1:4){ ni[[i]] <- sparseMatrixStats::rowCounts(x = vars_matrix, value = i-1) } # Careful the values are offset by one compared to vartrix conventions.
  ntot <- ni[[1]]+ni[[2]]+ni[[3]]+ni[[4]]
  for(i in 1:4){
    pi[[i]] <- ni[[i]]/ntot
    log2pi[[i]] <- pi[[i]] # Make a small work around to avoid log(0) -> if the probability is 0, p*log(p) is null, so log(p) can be set to 0 to get there.
    log2pi[[i]][pi[[i]]>0] <- log2(pi[[i]][pi[[i]]>0])
  }
  entropy <- -pi[[1]]*log2pi[[1]]-pi[[2]]*log2pi[[2]]-pi[[3]]*log2pi[[3]]-pi[[4]]*log2pi[[4]]
  names(entropy) <- rownames(vars_matrix)
  frac <- coverage/ncol(vars_matrix)
  minimal_entropy <- -frac*log2(frac)-(1-frac)*log2(1-frac)
  minimal_entropy[frac==0] <- 0
  excess_entropy <- entropy - minimal_entropy
  names(excess_entropy) <- rownames(vars_matrix)

  # Pack the entropy, coverage, and sorted variants in the genotype object:
  genotype$excess_entropy <- excess_entropy
  genotype$coverage <- coverage
  genotype$coverage_frac <- frac
  genotype@variants_by_coverage <- names(sort(-coverage))
  genotype@variants_by_information <- names(sort(-excess_entropy))
  genotype@informative_variants <- genotype@variants_by_information[1:n.variants]

  return(genotype)
}
