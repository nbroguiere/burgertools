#' Aggregate pseudo-bulk genotypes
#'
#' Starting from a Seurat object, for each of the identities given as the ident argument, compute a pseudo-bulk aggregated genotype. Store the results in a genotype object. If the ident is a list of cells, compute only one aggregated genotype for the cells listed. If it is the name of a metadata column, or a vector listing the identity of each individual cells, use the levels of this vector/column to define the groups of cells aggregated into individual pseudo-bulk genotypes.
#'
#' @param object Seurat object.
#' @param ident character(n). The name of a metadata column, or a vector listing the individual identities of cells, or a list of cells. Defines the groups of cells that are aggregated in pseudo-bulk genotypes. Default: Idents(object).
#' @param tolerance.pct numeric(1). The percentage of counts not matching the rest tolerated to keep a homozygous call.
#' @param assays character(3). The name of the assays in which the consensus matrix (expected to also contain variant metadata), reference allele count and alt allele count are stored in the seurat object. Default: c("VAR","REF","ALT").
#' @param slot character(1). The assay slot that should be used in the Seurat object. Default: "data".
#' @return Returns a genotype object with the aggregated pseudo-bulk genotypes and their metadata.
#' @keywords pseudo-bulk Pseudobulk genotype aggregation
#' @export
#' @examples
#' GenotypeObject <- PseudoBulkGenotypes(MySeuratObject)
PseudoBulkGenotypes <- function(object, ident=Idents(object), tolerance.pct=5, assays=c("VAR","REF","ALT"), slot="data"){
  # Clarify the identities to use
  if(!length(ident)){
    print("No cluster identity given, using default Idents(object)")
    idents <- Idents(object)
  }else{
    if(is.na(ident[1])){
      print("Using the levels of Idents(object) for pseudo-bulk genotype aggregation.")
      idents <- Idents(object)
    }else{
      if(length(ident)==1){
        print("Using metadata column",ident,"as identities for the cluster genotypes.")
        idents <- object@meta.data[,ident,drop=T]
      }else{
        if(length(intersect(ident,colnames(object)))==length(ident)){
          print("Using the list of cell barcodes given to aggregate a unique genotype.")
          object <- object[,ident]
          Idents(object) <- "GT"
        }else if(length(ident)==ncol(object)){
          if(sum(ident==Idents(object))==length(ident)){
            print("Using the levels of Idents(object) for pseudo-bulk genotype aggregation.")
            idents <- ident
          }else{
            print("Using the levels of the identities given to define the clusters aggregated into genotypes.")
            idents <- ident
          }
        }else{
          warning("Could not make sense of the ident given as input. Give a vector of cell names, a metadata column, a set of identities of the same length as ncol(object), or NA/NULL to use current Idents(object). Aborting.")
          return()
        }
      }
    }
  }
  # Names of the genotypes:
  gt.names <- as.character(unique(idents))

  # Pick up the variant metadata as is from the seurat object:
  meta <- object[[assays[1]]]@meta.features

  # Prepare a genotype sparse matrix to receive the data (vartrix conventions):
  m <- Matrix::Matrix(data = 0, nrow = nrow(object[[assays[1]]]), ncol = length(gt.names), dimnames = list(rownames(object[[assays[1]]]),gt.names))
  m <- methods::as(m,"dgCMatrix")

  # Go through the clusters, pick up the data, aggregate, compute coverage, frequency, and make a genotype call:
  tol <- tolerance.pct/100
  for(i in gt.names){
    ref <- Matrix::rowSums(GetAssayData(object, assay=assays[2], slot=slot)[,idents==i])
    alt <- Matrix::rowSums(GetAssayData(object, assay=assays[3], slot=slot)[,idents==i])
    cvg <- ref+alt
    freq <- alt/cvg + 1 # Offset by one to distinguish no data from 0 frequency
    freq[cvg==0] <- 0
    gt <- ref # Just use ref as a template for the shape
    gt[] <- 0 # No data unless overwritten
    gt[cvg>0] <- 3 # ref/alt unless overwritten
    gt[freq >=1 & freq<=1+tol] <- 1 # ref/ref
    gt[freq>=2-tol] <- 2 # alt/alt

    # fill in the matrices
    meta[,paste0("coverage_",i)] <- cvg
    m[,i] <- gt
    m <- Matrix::drop0(m)
  }

  # Package the results in a genotype object to return:
  return(CreateGenotypeObject(matrix = m, metadata = meta, variants = rownames(m), informative_variants = VariableFeatures(object, assay=assays[1])))
}
