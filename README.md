
# burgertools

The burgertools package contains various tools facilitating the analysis of samples with highly heterogeneous cell types. For example, it contains: 
- MAGIC imputation wrappers
- Signature scoring functions (best combined with imputed data)
- Cell type classifier based on markers (relying on signature scoring and MAGIC imputation to overcome dropouts for reliable classification based on even few markers)
- Cell type aware quality control (QC) plots.
- Quality control by cell type

## Installation

You can install the released version of burgertools from [github](https://github.com/nbroguiere/) with:

``` r
install.packages("githubinstall")
library(githubinstall)
githubinstall("burgertools")
```
or
``` r
install.packages("devtools")
library(devtools)
install_github("nbroguiere/burgertools")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(burgertools)
## basic example code
```

