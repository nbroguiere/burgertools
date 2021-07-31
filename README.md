
# burgertools

The burgertools package contains various tools facilitating the analysis of samples with highly heterogeneous cell types. It relies on Seurat for scRNA-seq data handling and MAGIC for imputation, and provides a user-friendly interface to classify cells based on signatures of just a few canonical markers, and make the data filtering by cell type. 
- Mitochondrial content determination and imputation with Seurat/MAGIC wrappers.
- Signature scoring functions (best combined with imputed data).
- Cell type classifier based on markers (relying on signature scoring and MAGIC imputation to overcome dropouts for reliable classification based on even few markers).
- Cell type aware quality control (QC) plots.
- Quality control by cell type.

## Installation

You can install the released version of burgertools from [github](https://github.com/nbroguiere/) with:

``` r
install.packages("devtools")
library(devtools)
install_github("nbroguiere/burgertools")
```
or
``` r
install.packages("githubinstall")
library(githubinstall)
githubinstall("burgertools")
```

## Example

Upcoming

``` r
library(burgertools)
## example
```
