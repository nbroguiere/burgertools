#' Annotate variants with VEP data.
#'
#' Add variant effect predictor (VEP) annotations to a genotype object, and renames the variants with a short human-intelligible description of the variant.
#'
#' If available in the VEP data frame, adds the fields: symbol, gene_ensembl_id, impact, protein_position, aminoacids, clinicalsign, summary, and- short_description to the genotype object metadata. Also renames the variants (typically named chr-pos-ref-alt after reading a vcf) with a short description prefix of the form symbol-protein_position-aminoacids-impact. Fields for which no data is available, or impact values less than min_impact, are skipped from the short description.
#'
#' VEP typically assigns several annotations for each variant in a vcf file. The one selected for annotating the genotype object is the one with the highest predicted impact. If several have the same impact, the ones with a symbol available are prioritized, and among those the ones with a protein position / amino acid change prediction. If there is still a draw between several annotations, the first one listed in the original VEP data is used.
#'
#'   refers to HGNC gene symbols.
#'
#' @param genotype A genotype object, typically from ReadVcf.
#' @param vep A data frame or tibble, typically from ReadVEP, containing the fields "Uploaded_variation", "Location", "Allele", "Gene", "SYMBOL", "IMPACT", "Consequence", "Protein_position", "Amino_acids", "CLIN_SIG"
#' @param min_impact character(1). Minimal impact included in the short descriptions of variants, among: "-"<"MODIFIER"<"LOW"<"MODERATE"<"HIGH". Default: "MODERATE"
#' @param avoid_underscores logical. Should underscores be replaced by dashes as a field separator. Needed when the variant data will be used with Seurat, which replaces underscore in feature names by dashes. Default: TRUE.
#' @return A genotype object with additional metadata fields (symbol gene_ensembl_id impact protein_position aminoacids clinicalsign consequence summary short_description) and annotated variant names (symbol-protein_position-aminoacids-impact-chr_pos-ref-alt).
#' @keywords annotate variants vep genotype
#' @examples
#' MyGenotype <- AnnotateWithVEP(MyGenotype,MyVepDataFrame,min.impact="HIGH",avoid_underscores=TRUE)
#' @export

AnnotateWithVEP <- function(genotype,vep,min_impact="MODERATE",avoid_underscores=TRUE){
  vep2 <- vep[vep$Uploaded_variation %in% genotype@metadata$VEP_ID,c("Uploaded_variation","Location","Allele","Gene","SYMBOL","IMPACT","Consequence","Protein_position","Amino_acids","CLIN_SIG")]
  annotations <- data.frame(VEP_ID=genotype@metadata$VEP_ID, symbol="-",gene_ensembl_id="-",impact="-",protein_position="-",aminoacids="-", consequence="-", clinicalsign="-",summary="-",short_description="-", row.names=genotype@metadata$UNIQUE_ID)
  annotations$impact <- factor(annotations$impact,levels=c("-","MODIFIER","LOW","MODERATE","HIGH"),ordered = T)

  # Create a score - 8 points per level of IMPACT, then one point for presence of a gene name, and one point for presence of the protein position.
  # Sort based on this, and then do a remove duplicates. Then build names vectorially.
  # Use this as a workaround to a for loop with if statements, that was running far too slow to be useful.
  vep2$priority <- 8*as.numeric(vep2$IMPACT) + as.numeric(vep2$SYMBOL!="-") + as.numeric(vep2$Protein_position!="-")
  vep2 <- vep2[order(-vep2$priority),]
  vep2 <- vep2[!duplicated(vep2$Uploaded_variation),]
  vep2 <- as.data.frame(vep2[,-ncol(vep2)])
  rownames(vep2) <- vep2$Uploaded_variation
  vep2$IMPACT2 <- vep2$IMPACT
  vep2$IMPACT2[vep2$IMPACT2<min_impact] <- "-"
  vep2$summary <- paste0("_",vep2$SYMBOL,"_",vep2$Protein_position,"_",vep2$Amino_acids,"_",vep2$IMPACT2)
  vep2$summary <- stringr::str_remove_all(vep2$summary,"_-")
  vep2$summary <- stringr::str_remove_all(vep2$summary,"^_")
  annotations_available <- intersect(annotations$VEP_ID,rownames(vep2))
  annotations$symbol[annotations$VEP_ID %in% annotations_available] <- vep2[annotations$VEP_ID[annotations$VEP_ID %in% annotations_available],"SYMBOL"]
  annotations$gene_ensembl_id[annotations$VEP_ID %in% annotations_available] <- vep2[annotations$VEP_ID[annotations$VEP_ID %in% annotations_available],"Gene"]
  annotations$impact[annotations$VEP_ID %in% annotations_available] <- vep2[annotations$VEP_ID[annotations$VEP_ID %in% annotations_available],"IMPACT"]
  annotations$protein_position[annotations$VEP_ID %in% annotations_available] <- vep2[annotations$VEP_ID[annotations$VEP_ID %in% annotations_available],"Protein_position"]
  annotations$aminoacids[annotations$VEP_ID %in% annotations_available] <- vep2[annotations$VEP_ID[annotations$VEP_ID %in% annotations_available],"Amino_acids"]
  annotations$clinicalsign[annotations$VEP_ID %in% annotations_available] <- vep2[annotations$VEP_ID[annotations$VEP_ID %in% annotations_available],"CLIN_SIG"]
  annotations$consequence[annotations$VEP_ID %in% annotations_available] <- vep2[annotations$VEP_ID[annotations$VEP_ID %in% annotations_available],"Consequence"]
  annotations$summary[annotations$VEP_ID %in% annotations_available] <- vep2[annotations$VEP_ID[annotations$VEP_ID %in% annotations_available],"summary"]
  annotations$short_description <- paste0(annotations$summary,"_",rownames(annotations))

  # Seurat doesn't accept underscores, so ends up in name conflicts later on. Replace underscores by dashes:
  if(avoid_underscores){
    rownames(annotations) <- stringr::str_replace_all(rownames(annotations),"_","-")
    annotations$summary <- stringr::str_replace_all(annotations$summary,"_","-")
    annotations$short_description <- stringr::str_replace_all(annotations$short_description,"_","-")
  }

  # Update the genotype object with the human understandable descriptive names:
  genotype@metadata <- cbind(annotations[,!colnames(annotations) %in% colnames(genotype@metadata)],genotype@metadata)
  rownames(genotype) <- annotations$short_description
  return(genotype)
}
