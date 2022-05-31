#' Read a vcf file
#'
#' This function reads a vcf file and returns the contents as a genotype object.
#'
#' This function reads a vcf file, containing the tab separated fields CHROM, POS, REF, ALT, ID, QUAL, INFO, FORMAT, and an undefined additional number of columns with genotypes (column titles are the genotype names, typically patient IDs), and returns a genotype object (genotypes as a sparse matrix in vartrix conventions, and associated metadata as a data frame).
#'
#' @param file The path to the vcf file, can be compressed (.vcf or .vcf.gz).
#' @return genotype object
#' @keywords genotypes variants vcf
#' @examples
#' genotypes <- ReadVcf("MyFolder/Myfile.vcf.gz")
#' @export

ReadVcf <- function(file=""){
  cat("Reading vcf file.\n")
  vcf <- readr::read_tsv(file = file, comment = "##",show_col_types = FALSE)
  colnames(vcf)[1] <- stringr::str_replace(colnames(vcf)[1],pattern = "#",replacement = "")

  cat("Converting to sparse column matrix in vartrix conventions.\n")
  GT_cols <- (which(colnames(vcf)=="FORMAT")+1):ncol(vcf)
  colnames(vcf)[GT_cols] <- paste0("Patient_",colnames(vcf)[GT_cols])
  GT <- vcf
  for(i in GT_cols){
    GT[,i] <- limma::strsplit2(dplyr::pull(vcf,var=i),":")[,1]
    cat(i-min(GT_cols)+1," out of ",length(GT_cols),"\n")
  }
  GT <- GT[,GT_cols]
  GT[GT=="./."] <- "0"
  GT[GT=="0/0"] <- "1"
  GT[GT=="1/1"] <- "2"
  GT[GT=="0/1"] <- "3"
  GT[GT=="1/0"] <- "3"
  GT <- as(apply(GT,2,as.numeric),"dgCMatrix")

  cat("Picking up the metadata separately, and removing the columns which contain no information.\n")
  meta <- vcf[,-GT_cols]
  tmp <- c()
  for(i in 1:ncol(meta)){
    if(length(unique(meta[,i,drop=T]))==1){
      tmp <- c(tmp,i)
    }
  }
  meta <- meta[,-tmp]

  cat("Generating unique row names of the form chr-position-ref-alt\n")
  UNIQUE_ID <- paste0(meta$CHROM,"-",format(meta$POS,scientific=FALSE,trim = TRUE),"-",meta$REF,"-",meta$ALT)
  rownames(GT) <- UNIQUE_ID
  colnames(GT) <- colnames(vcf[GT_cols])
  meta$UNIQUE_ID <- UNIQUE_ID

  cat("Generating (non-unique) IDs in VEP 'uploaded variant' conventions.\n")
  # Find the variant names in VEP "Uploaded variant" writing conventions
  # ("Uploaded variant" column in vep files - insertions do NOT have a REF allele, deletions do NOT have an ALT allele,
  # instead of showing one nucleotide being replaced. Additionally, this shifts the positions by one.
  meta$VEP_ID <- meta$ID # When there is a rs* ID present, it remains the VEP ID.
  no_id <- meta[meta$ID ==".",c("CHROM","POS","ID","REF","ALT","UNIQUE_ID")]
  keep_as_is_logical_indices <- (stringr::str_length(no_id$REF)==1 & stringr::str_length(no_id$ALT)==1) | no_id$REF=="-" | no_id$ALT=="-"
  keep_as_is <- no_id[keep_as_is_logical_indices,]
  keep_as_is$POS2 <- keep_as_is$POS
  keep_as_is$REF2 <- keep_as_is$REF
  keep_as_is$ALT2 <- keep_as_is$ALT
  need_shift <- no_id[!keep_as_is_logical_indices,]
  shift_left <- need_shift[stringr::str_length(need_shift$REF)==1,]
  shift_left$POS2 <- shift_left$POS+1
  shift_left$REF2 <- "-"
  shift_left$ALT2 <- substr(shift_left$ALT,2,stringr::str_length(shift_left$ALT))
  shift_right <- need_shift[stringr::str_length(need_shift$ALT)==1,]
  shift_right$POS2 <- shift_right$POS+1
  shift_right$REF2 <- substr(shift_right$REF,2,stringr::str_length(shift_right$REF))
  shift_right$ALT2 <- "-"
  no_id_vep_conventions <- rbind(shift_left,keep_as_is,shift_right)
  no_id_vep_conventions$VEP_ID <- paste0(no_id_vep_conventions$CHROM,"_",format(no_id_vep_conventions$POS2,scientific=FALSE,trim = TRUE),"_",no_id_vep_conventions$REF2,"/",no_id_vep_conventions$ALT2)
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$UNIQUE_ID
  meta[no_id_vep_conventions$UNIQUE_ID,"VEP_ID"] <- no_id_vep_conventions$VEP_ID

  cat("Returning genotype object\n")
  return(GenotypeObject(matrix=GT, metadata=meta, variants=meta$UNIQUE_ID))
}
