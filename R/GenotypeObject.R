#' Object storing genotypes and variant metadata
#'
#' The genotype objects are tailored to handle efficiently in R the contents of (mouse/human) vcf files, i.e. variant metadata, and several genotypes. They are also meant to be augmented with additional information such as annotations (e.g. vep) and single-cell analysis statistics.
#'
#' For efficient handling, the genotypes are stored as a sparse numeric matrix. When creating the object with ReadVcf, the vartrix conventions (from 10X genomics) are used, namely:
#' 0 for "no call" or in vcf "./."
#' 1 for "ref/ref" or in vcf "0/0"
#' 2 for "alt/alt" or in vcf "1/1"
#' 3 for "ref/alt" or in vcf "0/1"
#'
#' The objects can also be created manually typically with: GenotypeObject(matrix, metadata, variants).
#'
#' The genotype data is found in the matrix slot, and information concerning the variants is contained in the metadata slot. The list of variants currently in the object is found in the variants slot.
#'
#' Additional slots can be populated with additional data relevant to single cell analysis: most informative variants (informative_variants slot), and variants sorted by coverage (variants_by_coverage slot) or by entropy (variants_by_information).
#'
#' Standard generic methods can be used on the genotype object. In particular, the object can be subset with object[[i,j]] to obtain a genotype object with the corresponding subset of variants and genotypes (with full associated metadata, and restricted ranked variant lists). Object[i,j(,drop)] enables to access the genotype matrix in read and write, as well as rowSums, colSums, rownames, colnames, nrow, ncol, and dim. The $ operator enables to access directly to metadata columns, for both read and write.
#'
#' @slot matrix A dgCMatrix containing genotype data in vartrix conventions.
#' @slot metadata A data.frame containing variant metadata.
#' @slot variants A character vector listing the variants described in the object, in the same order as the rows of the matrix and metadata slots.
#' @slot vartrix A list of vartrix matrices
#' @slot variants_by_coverage A character vector listing the variants sorted from max coverage (number of cells in which there is data) to min.
#' @slot variants_by_information A character vector listing the variants sorted from max information (excess entropy in single cell data) to min.
#' @slot informative_variants A character vector listing the most informative variants, used for downstream clustering analysis.
#' @keywords genotype vcf vep variants
#' @export
#' @examples
#' genotypes <- ReadVcf("Myfolder/MyVcf.vcf")
#' genotypes # Preview the contents of the genotype object
#' genotypes[1:3,1:6] # Check the genotypes matrix
#' genotypes[1:2,1:5] <- 2 # Modify the genotypes matrix
#' head(genotypes@metadata) # Check the genotypes metadata
#' genotypes$CHROM[1:5] <- "chrZ" # Modify a metadata column
#' genotypes$CHROM[1:10] # Retrieve elements of a metadata column
#' genotypes$CHROM <- NULL # Delete a metadata column
#' colnames(genotypes)[1:3] <- c("a","b","c") # Rename some genotypes
#' colnames(genotypes) # Access genotype names
#' rownames(genotypes)[1:2] <- c("a","b") # Rename some variants
#' nrow(genotypes) # Get the number of variants
#' ncol(genotypes) # Get the number of genotypes
#' dim(genotypes) # Get both
#' genotypes@informative_variants <- genotypes@variants[1:20] # Set the informative variants list
#' genotypes[[1:10,1:5]] # Subset the genotype object
#' genotypes[[rownames(genotypes)[1:5],colnames(genotypes)[1:5]]] # Subset the object by variant and patient name (returns a genotype object)
#' genotypes[rownames(genotypes)[1:5],colnames(genotypes)[1:5]] # Access the genotype matrix by variant and patient name
GenotypeObject <- setClass("genotype", slots=list(matrix="dgCMatrix",
                                                        metadata="data.frame",
                                                        variants="character",
                                                        variants_by_coverage="character",
                                                        variants_by_information="character",
                                                        informative_variants="character",
                                                        vartrix="list"))

#' @describeIn GenotypeObject Show a summary of the contents of a genotype object.
#' @export
setMethod("show","genotype",
          function(object){
            cat("A genotype S4 object containing ",length(object@variants)," variants (@variant slot) for ",ncol(object@matrix)," unique genotypes (@matrix slot) and ",ncol(object@metadata)," variant annotations (@metadata slot).\n\n")
            cat("Head(3) of the genotype matrix:\n\n")
            print(head(object@matrix,3))
            cat("\n\nHead(3) of the metadata data frame:\n\n")
            print(head(object@metadata,3))
            tmp <- ""
            if(length(object@variants_by_coverage)) tmp <- paste0(tmp,"@variants_by_coverage (",length(object@variants_by_coverage),")\n")
            if(length(object@variants_by_information)) tmp <- paste0(tmp,"@variants_by_information (",length(object@variants_by_information),")\n")
            if(length(object@informative_variants)) tmp <- paste0(tmp,"@informative_variants (",length(object@informative_variants),")\n")
            if(length(object@vartrix)) tmp <- paste0(tmp,"@vartrix (",paste0(names(object@vartrix),collapse = ", "),")\n")
            if(tmp!=""){
              cat(paste0("\n\nOther slots populated:\n",tmp))
            }
          }
)

#' @describeIn GenotypeObject Access matrix values in a genotype object.
#' @export
setMethod("[","genotype",
          function(x, i, j, ..., drop=F){
            if(missing(i)) i=T
            if(missing(j)) j=T
            return(x@matrix[i,j, ..., drop])
          }
)

#' @describeIn GenotypeObject Assign genotype matrix values in a genotype object
#' @export
setMethod("[<-","genotype",
          function(x, i, j, ..., value){
            if(missing(i)) i=T
            if(missing(j)) j=T
            x@matrix[i,j] <- value
            return(x)
          }
)

#' @describeIn GenotypeObject Access metadata columns in a genotype object
#' @export
setMethod("$","genotype",
          function(x,name){
            tmp <- x@metadata[,name,drop=T]
            names(tmp) <- x@variants
            return(tmp)
          }
)

#' @describeIn GenotypeObject Subset a genotype object
#' @export
setMethod("[[","genotype",
          function(x, i, j, ...){
            if(missing(i)) i=T
            if(missing(j)) j=T
            y <- GenotypeObject(
              matrix=x@matrix[i,j,drop=F],
              metadata=x@metadata[i,,drop=F],
              variants=stats::setNames(x@variants,x@variants)[i],
              vartrix=x@vartrix
            )
            if(length(x@variants_by_coverage)) y@variants_by_coverage=intersect(x@variants_by_coverage,y@variants)
            if(length(x@variants_by_information)) y@variants_by_information=intersect(x@variants_by_information,y@variants)
            if(length(x@informative_variants)) y@informative_variants=intersect(x@informative_variants,y@variants)
            if(length(y@vartrix)){
              for(k in 1:length(y@vartrix)){
                y@vartrix[[k]] <- y@vartrix[[k]][y@variants,]
              }
            }
            return(y)
          }
)

#' @describeIn GenotypeObject Assign values to genotype object metadata columns
#' @export
setMethod("$<-","genotype",
          function(x,name,value){
            if(!name %in% colnames(x@metadata) & !is.null(value)){
              x@metadata <- cbind(x@metadata,value)
              colnames(x@metadata)[ncol(x@metadata)] <- name
            }else{
              x@metadata[,name] <- value
            }
            return(x)
          }
)

.DollarNames.genotype <- function(x,pattern=""){
  grep(pattern, colnames(x@metadata), value=TRUE)
}

#' @describeIn GenotypeObject Retrieve column names (i.e. genotype names) from a genotype object.
#' @export
setMethod("colnames","genotype",
          function(x){
            return(colnames(x@matrix))
          }
)

#' @describeIn GenotypeObject Assign column names (i.e. genotype names) to a genotype object.
#' @export
setMethod("colnames<-","genotype",
          function(x,value){
            colnames(x@matrix) <- value
            return(x)
          }
)

#' @describeIn GenotypeObject Retrieve row names (i.e. variant names) from a genotype object.
#' @export
setMethod("rownames","genotype",
          function(x){
            return(rownames(x@matrix))
          }
)

#' @describeIn GenotypeObject Assign row names (i.e. variant names) from a genotype object.
#' @export
setMethod("rownames<-","genotype",
          function(x,value){
            rownames(x@matrix) <- value
            rownames(x@metadata) <- value
            x@variants <- value
            if(length(x@vartrix)){
              for(i in 1:length(x@vartrix)){
                rownames(x@vartrix[[i]]) <- value
              }
            }
            return(x)
          }
)

#' @describeIn GenotypeObject Retrieve the number of rows (i.e. number of variants) from a genotype object.
#' @export
setMethod("nrow","genotype",
          function(x){
            return(nrow(x@matrix))
          }
)

#' @describeIn GenotypeObject Retrieve the number of columns (i.e. number of genotypes) from a genotype object.
#' @export
setMethod("ncol","genotype",
          function(x){
            return(ncol(x@matrix))
          }
)

#' @describeIn GenotypeObject Retrieve the dimension, i.e. number of variants and samples, from a genotype object.
#' @export
setMethod("dim","genotype",
          function(x){
            return(dim(x@matrix))
          }
)

#' @describeIn GenotypeObject Compute the number of genotypes which cover a given variant, from a genotype object.
#' @export
setMethod("rowSums","genotype",
          function(x){
            return(rowSums(x@matrix>0))
          }
)

#' @describeIn GenotypeObject Compute the number of variants covered in each genotype, from a genotype object.
#' @export
setMethod("colSums","genotype",
          function(x){
            return(colSums(x@matrix>0))
          }
)
