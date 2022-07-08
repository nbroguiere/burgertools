#' Load vartrix data
#'
#' Read the vartrix REF and ALT matrices, and computes frequency (FREQ) and consensus (VAR) matrices. Adds them to a genotype object already containing the reference vcf data.
#'
#' The FREQ sparse matrix is offset by 1 in order to efficiently distinguish zero-frequency (value 1) from missing data (value 0, not stored in the sparse matrix).
#'
#' The consensus matrix makes a call for no data (value 0, not stored in the sparse matrix), ref/ref (value 1), alt/alt (value 2) or alt/ref (value 3) genotypes. Since in scRNA-seq data, there is likely some extracellular RNA from dead cells or debris slightly contaminating the reads of other cells, a tolerance can be set above 0. For example, at a tolerance of 5%, up to 1 read in 20 can differ from the rest without making the call switch to heterozygous. The behavior of the original vartrix consensus matrix calculation corresponds to a tolerance of 0.
#'
#' Equivalence between vartrix, human, and vcf genotype naming conventions:
#'
#' "0" for "no call", "./."
#'
#' "1" for "ref/ref", "0/0"
#'
#' "2" for "alt/alt", "1/1"
#'
#' "3" for "ref/alt", "0/1"
#'
#' @param genotype A genotype object, which contains the reference vcf file used to compute vartrix matrices.
#' @param ref character(1), the file name for the reference allele count matrix.
#' @param alt character(1), the file name for the alternative allele count matrix.
#' @param barcodes character(n), the file name for the cell barcodes, or the list of cell barcodes for which vartrix was run, used as column names on the matrices. c() for no names. Default: c().
#' @param tolerance.percent numeric(1), the percentage of alt or ref counts not in agreement with the rest on a unique (cell,variant) tolerated without affecting the call. Default: 5%.
#' @param strip.suffix logical(1), should the "-X" suffix be striped from the barcodes if constant. Default: TRUE.
#' @return Returns the genotype object with a populated vartrix slot, containing a list of sparse matrices (dgCmatrix), with names REF, ALT, FREQ and VAR (=consensus). Variants are the rows and cells/barcodes the columns.
#' @keywords vartrix variants matrix ref alt consensus freq frequency genotype
#' @examples
#' MyGenotypes <- ReadVartrix(MyGenotypes, "vartrix_ref_matrix.mtx.gz", "vartrix_alt_matrix.mtx.gz", barcodes=MyCellBarcodes, tolerance.percent=2)
#' @export

ReadVartrix <- function(genotype, ref, alt, barcodes=c(), tolerance.percent=5, strip.suffix=TRUE){
  if(length(barcodes)==0){
    warning("No cell barcodes given to ReadVartrix. Beware that the column names for vartrix matrices will not be set.\n")
  }
  if(length(barcodes)==1){
    cat("Reading barcodes from file.\n")
    barcodes <- readr::read_tsv(barcodes,col_types = "c",show_col_types = F,col_names = F)[,1,drop=T]
    if(strip.suffix){
      tmp <- limma::strsplit2(barcodes,"-")
      if(length(unique(tmp[,2,drop=T]))==1){
        barcodes <- tmp[,1,drop=T]
      }
    }
  }
  if(!length(genotype@variants)){
    warning("The genotype does not contain variants (genotype@variants not populated), the row names for the vartrix matrices will not be set.
            You can set rownames with 'rownames(genotype_object)<-' but it is recommended to rather read the reference vcf file into the genotype object first.")
  }

  cat("Reading the reference allele count matrix\n")
  vartrix_ref <- Matrix::drop0(Matrix::readMM(ref))
  if(dim(vartrix_ref)[1]!=length(genotype@variants)){
    warning(paste0("The number of lines in the vartrix reference matrix (",dim(vartrix_ref)[1],") does not match the number of variants in the genotype object (",length(genotype@variants),"). Aborting."))
    return(genotype)
  }
  if(dim(vartrix_ref)[2]!=length(barcodes)){
    warning(paste0("The number of columns in the vartrix reference matrix (",dim(vartrix_ref)[2],") does not match the number of barcodes (",length(barcodes),"). Aborting."))
    return(genotype)
  }
  if(length(barcodes))          colnames(vartrix_ref) <- barcodes
  if(length(genotype@variants)) rownames(vartrix_ref) <- genotype@variants

  cat("Reading the alternate allele count matrix\n")
  vartrix_alt <- Matrix::drop0(Matrix::readMM(alt))
  if(length(barcodes))          colnames(vartrix_alt) <- barcodes
  if(length(genotype@variants)) rownames(vartrix_alt) <- genotype@variants

  cat("Computing the alt frequency matrix (offset by one to distinguish 0 frequency from no data in the sparse matrix) \n")
  tmp <- vartrix_ref+vartrix_alt
  tmp@x <- 1/tmp@x
  vartrix_freq_matrix <- 0*vartrix_ref+vartrix_alt*tmp # In this sparse matrix, explicit zeros are actual zero frequency of alt, missing values mean no data.
  vartrix_freq_matrix@x <- 1 + vartrix_freq_matrix@x # Offset all the data by one so that zeros from no data (not shifted) can be distinguished from zero frequency (shifted to one)

  # Compute the consensus matrix, use the vartrix conventions:
  # "0" for "no call", "./."
  # "1" for "ref/ref", "0/0"
  # "2" for "alt/alt", "1/1"
  # "3" for "ref/alt", "0/1"
  # Noise tolerant version, necessary to avoid that a small amount of ambiant RNA or debris, especially in very high expressed genes (e.g. mitochondrial), throws off the calculation:
  cat("Computing the consensus matrix, with an error tolerance of ",tolerance.percent,"%\n")
  vartrix_consensus_matrix <- vartrix_freq_matrix
  vartrix_consensus_matrix@x[vartrix_freq_matrix@x <= 1+tolerance.percent/100] <- 1
  vartrix_consensus_matrix@x[vartrix_freq_matrix@x >= 2-tolerance.percent/100] <- 2
  vartrix_consensus_matrix@x[vartrix_freq_matrix@x > 1+tolerance.percent/100 & vartrix_freq_matrix@x < 2-tolerance.percent/100] <- 3

  # Of note the algorithm above is universal, but if the tolerance is zero (as is the case for vartrix consensus), I could have used the simpler and very elegant code below.
  # Runs so fast, let's actually run it anyway just to compare for curiosity the proportion of the calls that were affected by the tolerance.
  # 0.45% in a test run with 8 multiplexed patients and a lot of contaminating RNA, for a tolerance of 7%.
  vartrix_consensus_matrix_traditional <- 1*(vartrix_ref>0)+2*(vartrix_alt>0)
  cat("Percentage of calls corrected because of the tolerance: ", 100*sum(abs(vartrix_consensus_matrix-vartrix_consensus_matrix_traditional))/sum(abs(vartrix_consensus_matrix)),"%\n")

  genotype@vartrix <- list(REF=vartrix_ref, ALT=vartrix_alt, FREQ=vartrix_freq_matrix, VAR=vartrix_consensus_matrix)
  return(genotype)
}

