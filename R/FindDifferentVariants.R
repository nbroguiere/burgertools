#' Find differentially frequent variants between two populations
#'
#' Starting from a Seurat object, and two groups of cells (given as their ident name, or as a vector of cell barcodes. If second group is omitted, assumed to be all other cells), computes the differential frequency of occurance of variants with a chi-square frequency test. Both p-values and Bonferroni-corrected adjusted p-values are returned, as well as the difference of frequency of alternate alleles between the two populations.
#'
#' @param object Seurat object.
#' @param ident.1 character(n). The name of an identity or a list of cell barcodes or column positions in the Seurat object, defining the first cell population.
#' @param ident.2 character(n). The name of an identity or a list of cell barcodes or column positions in the Seurat object, defining the second cell population, or NULL to use all the other cells. Default: NULL.
#' @param min.cells numeric(1). The minimal number of cells that should be have data in both populations for differential testing to be performed. Default: 15.
#' @param sort.by character(1). "pval" to sort the result dataframe by p-value, or "deltafreq" to sort by absolute difference in frequency. Any other value will result in no sorting, conserving the feature order in the Seurat object. Default: pval.
#' @param assays character(2). The name of the assays in which the reference allele count and alt allele count are stored in the seurat object. Default: c("REF","ALT").
#' @param layer character(1). The assay layer that should be used in the Seurat object. Default: "data".
#' @param show.progress numeric(1). Number of tests performed between updates on the advance. FALSE or 0 for no progress reports. Default: 500.
#' @param pval.lower.limit numeric(1). p-vals lower than this limit will be truncated to this limit. Useful to avoid log(0) and set an axis limit in volcano-like plots of the results.
#' @return Returns a data frame with columns "n1" (number of cells with variant data in population 1) and "n2" (same in population 2), the frequency of ref and alt allele in both populations "freq1" and "freq2", the difference in frequency "deltafreq", and the results of the chi-square test for identitical frequency with or without Bonferroni correction "pval" and "adj.pval".
#' @keywords differentially frequent variants statistical testing
#' @export
#' @examples
#' DifferentialTestingresults <- FindDifferentVariants(MySeuratObject, "FirstIdentName", "OtherIdentName")
FindDifferentVariants <- function(object, ident.1, ident.2=NULL, min.cells=15, sort.by="pval", assay=c("REF","ALT"), layer="data", show.progress=500, pval.lower.limit=1e-100){
  # Clarify the identities to use
  if(length(ident.1)==1){
    ident.1 <- colnames(object)[Idents(object)==ident.1]
  }
  if(length(ident.2)==1){
    ident.2 <- colnames(object)[Idents(object)==ident.2]
  }else if(is.null(ident.2)){
    ident.2 <- setdiff(colnames(object),ident.1)
  }

  if(!length(ident.1)){
    warning("ident.1 is empty. Aborting.")
    return(object)
  }
  if(!length(ident.2)){
    warning("ident.2 is empty. Aborting.")
    return(object)
  }
  # Pick up the data filtered above min.cell threshold, and aggregate the counts of ref and alt:
  variants_ident_1_ref <- GetAssayData(object, assay=assay[1], layer=layer)[,ident.1]
  variants_ident_1_alt <- GetAssayData(object, assay=assay[2], layer=layer)[,ident.1]
  variants_ident_2_ref <- GetAssayData(object, assay=assay[1], layer=layer)[,ident.2]
  variants_ident_2_alt <- GetAssayData(object, assay=assay[2], layer=layer)[,ident.2]
  tmp1 <- Matrix::rowSums((variants_ident_1_ref+variants_ident_1_alt)>0)
  tmp2 <- Matrix::rowSums((variants_ident_2_ref+variants_ident_2_alt)>0)
  keep <- tmp1>min.cells & tmp2>min.cells
  if(!sum(keep)){
    stop("No variants are present in the minimum requested number of cells in both populations, aborting.")
  }
  variant.names <- rownames(variants_ident_1_ref)[keep]
  variants_ident_1_ref <- Matrix::rowSums(variants_ident_1_ref[keep,])
  variants_ident_1_alt <- Matrix::rowSums(variants_ident_1_alt[keep,])
  variants_ident_2_ref <- Matrix::rowSums(variants_ident_2_ref[keep,])
  variants_ident_2_alt <- Matrix::rowSums(variants_ident_2_alt[keep,])

  # Prepare the dataframe that will collect the differential testing results
  difftest <- data.frame(row.names = variant.names)
  difftest$n1 <- NA
  difftest$n2 <- NA
  difftest$freq1 <- NA
  difftest$freq2 <- NA
  difftest$deltafreq <- NA
  difftest$pval <- NA
  difftest$adj.pval <- NA

  # Loop through the variants and perform chisquare tests
  if(show.progress){
    cat("Progress:\n")
  }
  for(i in 1:nrow(difftest)){
    xref <- variants_ident_1_ref[i]
    xalt <- variants_ident_1_alt[i]
    yref <- variants_ident_2_ref[i]
    yalt <- variants_ident_2_alt[i]
    contingency_table <- matrix(c(xref,xalt,yref,yalt), nrow = 2, ncol = 2, byrow = T, dimnames = list(c("x","y"),c("REF","ALT")))
    Xsq <- stats::chisq.test(contingency_table, correct = F)

    difftest$n1[i] <- xref+xalt
    difftest$n2[i] <- yref+yalt
    difftest$freq1[i] <- xalt/(xref+xalt)
    difftest$freq2[i] <- yalt/(yref+yalt)
    difftest$pval[i] <- Xsq$p.value

    if(show.progress){
      if(!i%%show.progress){
        print(paste0(i," out of ",nrow(difftest)))
      }
    }
  }
  difftest$deltafreq <- difftest$freq1 - difftest$freq2
  difftest$adj.pval  <- pmin(1,difftest$pval * nrow(difftest))

  difftest$pval     <- pmax(pval.lower.limit,difftest$pval)
  difftest$adj.pval <- pmax(pval.lower.limit,difftest$adj.pval)

  if(sort.by[1]=="deltafreq"){
    difftest <- difftest[order(-abs(difftest$deltafreq)),] # sort by difference in frequency
  }else if(sort.by[1]=="pval"){
    difftest <- difftest[order(difftest$pval),] # sort by significance
  }

  return(difftest)
}
