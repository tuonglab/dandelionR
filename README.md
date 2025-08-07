[![codecov](https://codecov.io/gh/tuonglab/dandelionR/graph/badge.svg?token=dd1bBTW48K)](https://codecov.io/gh/tuonglab/dandelionR)
[![R-CMD-check](https://github.com/tuonglab/dandelionR/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/tuonglab/dandelionR/actions/workflows/R-CMD-check.yml)
[![vignette](https://github.com/tuonglab/dandelionR/actions/workflows/vignette.yml/badge.svg)](https://github.com/tuonglab/dandelionR/actions/workflows/vignette.yml)

# dandelionR

Welcome to `dandelionR`!

`dandelionR` is an R package for performing single-cell immune repertoire trajectory analysis, based on the original python implementation in [`dandelion`](https://www.github.com/zktuong/dandelion).

It provides all the necessary tools to interface with [`scRepertoire`](https://github.com/ncborcherding/scRepertoire) and a custom implementation of absorbing markov chain for pseudotime inference, inspired based on the [palantir](https://github.com/dpeerlab/Palantir) python package.

## Installation from Bioconductor

You can install `dandelionR` from Bioconductor with:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# If installing from Bioc devel, do this first:
# BiocManager::install(version="devel")

BiocManager::install("dandelionR")
```

## Installation from GitHub

You can install `dandelionR` from GitHub with:

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github('tuonglab/dandelionR', dependencies = TRUE)
```

## Quick Start

```R
library(dandelionR)
```

This is a work in progress, so please feel free to open an issue if you encounter any problems or have any suggestions for improvement.

## Citation

If you use `dandelionR` in your work, you can cite our manuscript published in the [Computational and Structural Biotechnology Journal](https://authors.elsevier.com/sd/article/S2001-0370(25)00263-6):

```
Yu, J., Xu, X., Borcherding, N., Tuong, Z.K. dandelionR: Single-cell immune repertoire trajectory analysis in R. Comput. Struct. Biotechnol. J. (2025). doi:10.1016/j.csbj.2025.06.047.
```

Please also cite the original `dandelion` paper:

```
Suo, C. et al. Dandelion uses the single-cell adaptive immune receptor repertoire to explore lymphocyte developmental origins. Nat. Biotechnol. 42, 40-51 (2024). doi:10.1038/s41587-023-01734-7
```
