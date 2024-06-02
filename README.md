# Assessing Immune Factors in Maternal Milk and Paired Infant Plasma Antibody Binding to Human Rhinoviruses





The code provided here will process the peptide array data available under 10.5281/zenodo.8269759
and return plasma IgG epitope calls and breast milk IgA epitope calls which are 
the same as the epitope calls provided as supplemental data with our manuscript.

## Step-by-step directions

1. Clone this repo.
2. Download the input data (IgG.lseq.tsv and IgA.lseq.tsv) from Zenodo at https://zenodo.org/records/8269759 and 
add to this repo directory.
3. Run `epitope_calls_analysis.R`.

## Dependencies

This code was tested with `R-4.0.5` and `R 4.1.0` and requires the 
`dplyr`, `Matrix`, `data.table`, `openxlsx`, `tidyr`, and `purrr` packages and 
their dependencies. 

An example `R-4.0.5` session info is provided below:

```
 version  R version 4.0.5 (2021-03-31)
 os       macOS 13.4.1
 system   x86_64, darwin17.0

Packages ───────────────────────────────────────────────────────────────────────
 package         version
 circlize        0.4.15
 cli             3.6.1
 colorspace      2.0-3
 data.table      1.14.4
 dichromat       2.0-0.1
 digest          0.6.29
 dplyr           1.1.2
 ellipsis        0.3.2
 evaluate        0.17
 fansi           1.0.3
 fastmap         1.1.0
 forcats         0.5.1
 generics        0.1.3
 GlobalOptions   0.1.2
 glue            1.6.2
 htmltools       0.5.5
 httr            1.4.4
 knitr           1.40
 lattice         0.20-45
 lifecycle       1.0.3
 magrittr        2.0.3
 mapproj         1.2.8
 maps            3.4.0
 Matrix          1.5-3
 openxlsx        4.2.5.1
 pals            1.7
 pillar          1.9.0
 pkgconfig       2.0.3
 purrr           0.3.4
 qs              0.25.3
 R6              2.5.1
 RApiSerialize   0.1.2
 Rcpp            1.0.8
 RcppParallel    5.1.5
 rlang           1.1.1
 rmarkdown       2.17
 rstudioapi      0.14
 sessioninfo     1.2.2
 shape           1.4.6
 stringfish      0.15.5
 stringi         1.7.6
 stringr         1.4.1
 tibble          3.2.1
 tidyr           1.2.0
 tidyselect      1.2.0
 utf8            1.2.2
 vctrs           0.6.3
 withr           2.5.0
 xfun            0.39
 yaml            2.3.5
 zip             2.2.2
```
