#' Load vartrix data
#'
#' Read the vartrix ref and alt matrices, and computes a frequency and a consensus matrix. Returns all of them as a list of sparse matrices.
#'
#' The frequency sparse matrix is offset by 1 in order to distinguish zero-frequency (value 1) from missing data (value 0, not stored in the sparse matrix).
#'
#' The consensus matrix makes a call for no data (value 0, not stored in the sparse matrix), ref/ref (value 1), alt/alt (value 2) or alt/ref (value 3) genotypes. Since in scRNA-seq data, there is likely some extracellular RNA from dead cells or debris slightly contaminating the reads of other cells, a tolerance can be set above 0. For example, at a tolerance of 5%, up to 1 read in 20 can differ from the rest without making the call switch to heterozygous. The behavior of the vartrix consensus matrix calculation corresponds to a tolerance of 0.
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
#' @param ref_matrix character(1), the file name for the reference allele count matrix.
#' @param alt_matrix character(1), the file name for the alternate allele count matrix.
#' @param cell_barcodes character vector, the list of cell barcodes for which vartrix was run, used as column names on the matrices.
#' @param variants character vector, the list of variants for which vartrix was run, used as row names on the matrices.
#' @param tolerance_pct numeric(1), the percentage of alt or ref counts not in agreement with the rest on a unique (cell,variant) tolerated without affecting the call. Default: 5%.
#' @return A list of sparse matrices (dgCmatrix), with names ref, alt, freq, and consensus, variants as rows and cells/barcodes as columns.
#' @keywords vartrix variants matrix ref alt consensus freq frequency
#' @examples
#' vartrix_matrices <- ReadVartrix("MyVartrixRef.mtx.gz","MyVartrixAlt.mtx.gz",MyCellBarcodes,MyVariants)
#' @export

ReadVartrix <- function(ref_matrix,alt_matrix,cell_barcodes,variants,tolerance_pct=5){
  cat("Reading the reference allele count matrix\n")
  vartrix_ref_matrix <- Matrix::drop0(readMM(ref_matrix))
  colnames(vartrix_ref_matrix) <- cell_barcodes
  rownames(vartrix_ref_matrix) <- variants

  cat("Reading the alternative allele count matrix\n")
  vartrix_alt_matrix <- Matrix::drop0(readMM(alt_matrix))
  colnames(vartrix_alt_matrix) <- cell_barcodes
  rownames(vartrix_alt_matrix) <- variants

  cat("Computing the alt frequency matrix (offset by one to distinguish 0 frequency from no data in the sparse matrix) \n")
  tmp <- vartrix_ref_matrix+vartrix_alt_matrix
  tmp@x <- 1/tmp@x
  vartrix_freq_matrix <- 0*vartrix_ref_matrix+vartrix_alt_matrix*tmp # In this sparse matrix, explicit zeros are actual zero frequency of alt, missing values mean no data.
  vartrix_freq_matrix@x <- 1 + vartrix_freq_matrix@x # Offset all the data by one so that zeros from no data (not shifted) can be distinguished from zero frequency (shifted to one)
  #hist.data = hist((vartrix_freq_matrix@x-1)*100, 500, plot=F); hist.data$counts = log10(hist.data$counts); plot(hist.data) # histogram of the frequency matrix, helps to choose the tolerance.

  # Compute the consensus matrix, use the vartrix conventions:
  # "0" for "no call", "./."
  # "1" for "ref/ref", "0/0"
  # "2" for "alt/alt", "1/1"
  # "3" for "ref/alt", "0/1"
  # Noise tolerant version, necessary to avoid that a small amount of ambiant RNA or debris, especially in very high expressed genes (e.g. mitochondrial), throws off the calculation:
  cat("Computing the consensus matrix, with an error tolerance of ",tolerance_pct,"%\n")
  vartrix_consensus_matrix <- vartrix_freq_matrix
  vartrix_consensus_matrix@x[vartrix_freq_matrix@x <= 1+tolerance_pct/100] <- 1
  vartrix_consensus_matrix@x[vartrix_freq_matrix@x >= 2-tolerance_pct/100] <- 2
  vartrix_consensus_matrix@x[vartrix_freq_matrix@x > 1+tolerance_pct/100 & vartrix_freq_matrix@x < 2-tolerance_pct/100] <- 3

  # Of note the algorithm above is universal, but if the tolerance is zero (as is the case for vartrix consensus), I could have used the simpler and very elegant code below.
  # Runs so fast, let's actually run it anyway just to compare for curiosity the proportion of the calls that were affected by the tolerance.
  # 0.45% in a test run with 8 multiplexed patients and a lot of contaminating RNA, for a tolerance of 7%.
  vartrix_consensus_matrix_traditional <- 1*(vartrix_ref_matrix>0)+2*(vartrix_alt_matrix>0)
  cat("Percentage of calls affected by the tolerance: ", 100*sum(abs(vartrix_consensus_matrix-vartrix_consensus_matrix_traditional))/sum(abs(vartrix_consensus_matrix)),"%\n")

  return(list(ref=vartrix_ref_matrix, alt=vartrix_alt_matrix, freq=vartrix_freq_matrix, consensus=vartrix_consensus_matrix))
}

