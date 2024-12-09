[![codecov](https://codecov.io/gh/tuonglab/dandelionR/graph/badge.svg?token=dd1bBTW48K)](https://codecov.io/gh/tuonglab/dandelionR)

# dandelionR

Welcome to `dandelionR`!

`dandelionR` is an R package for performing single-cell immune repertoire trajectory analysis, based on the original python implementation in [`dandelion`](https://www.github.com/zktuong/dandelion).

It provides all the necessary tools to interface with [`scRepertoire`](https://github.com/ncborcherding/scRepertoire) and a custom implementation of absorbing markov chain for pseudotime inference, inspired based on the [palantir](https://github.com/dpeerlab/Palantir) python package.

## Installation

You can install `dandelionR` from GitHub with:

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
devtools::install_github('tuonglab/dandelionR', dependencies = TRUE)
```

## Quick Start

```R
library(dandelionR)
```

This is a work in progress, so please feel free to open an issue if you encounter any problems or have any suggestions for improvement.
