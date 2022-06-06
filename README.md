
# burgertools

The burgertools package contains various tools facilitating the analysis of tumor or other samples with highly heterogeneous cell types. It is built as an extension of Seurat. The package implements tools that facilitate interactive visualization of cell heterogeneity in FACS-like scatter plots (based on cell-type markers, scored as signatures on imputed data, or mutations/SNPs/variants, or cell hashing), and definition of gates in natural language to separate subcategories (e.g. defining cell types, interactively performing genetic demultiplexing, separating cancer cells from somatic cells based on variants. This is very useful when facing complex samples in order to perform QC by cell type, separate multiplexed samples (with both genetic and hashing-based combined multiplexing), or separate cell categories based on a limited number of canonical markers (e.g. CD3+ T cells, EPCAM+ epithelial cells, CD79B+ B cells, COL1A1+ fibroblasts), all at the single cell level and without clustering. 

In particular, the package includes:
- Wrappers for imputation with MAGIC, for Seurat objects, with sensible default parameters and result storage.
- Single-cell signature scoring functions (best combined with imputed data to overcome dropouts).
- Visualization of of signatures in a FACS-like manner in interactive scatter plots, QC plots or Feature plots, used to interactively define thresholds/gates. 
- Cell type classifier based on signatures, either automatic or gate-based and user-defined.
- Cell type aware quality control, essential for samples with very high cellular heterogeneity (e.g. T cells and cancer cells).
- Import and Export of filtered data and metadata in a 10X-like, open format (tsv/mtx), directly to and from Seurat objects. 
- Import of variant data directly from vcf files, to dedicated genotype objects easy to handle in R, and used to aggregate variant information.
- Import of variant annotations from TCGA and VEP, and integration with vcf data in genotype objects.
- Import of single-cell level variant data, obtained directly on scRNA-seq with vartrix, into Seurat objects. 
- Integration of genotype and Seurat objects to find informative variants, filter variants based on single-cell statistics as well as annotations, and plot them at the single cell level. 
- Fast dimensionality reductions based on variant data, for interactive exploration of genetic/chromosome silencing diversity. 
- Demultiplexing based on variants (blind or with reference WGS/WES patient-specific data), or analysis of tumor heterogeneity based on variants, with an interactive visual and gate-based strategy. 
- Demultiplexing based on hashtags, with a visual gate-based approach. 

## Installation

Install from [github](https://github.com/nbroguiere/burgertools) with:

``` r
install.packages("devtools")
devtools::install_github("nbroguiere/burgertools")
```
