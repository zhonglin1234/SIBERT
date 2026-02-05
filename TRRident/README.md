# TRRident <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/yourusername/TRRident/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/TRRident/actions)
<!-- badges: end -->

## T Cell Entry Region Identification in Spatial Transcriptomics

TRRident is a computational framework for identifying T cell entry regions (TRR) in spatial transcriptomics data. The method uses Ising models with impurities to detect regions where T cells enter the tumor core.

## Installation

You can install TRRident from GitHub:

```r
# install.packages("devtools")
devtools::install_github("yourusername/TRRident")
```

### Dependencies

TRRident requires:
- R (>= 4.0.0)
- Seurat (>= 4.0.0)
- Rcpp and RcppArmadillo
- BayesLogit
- mclust
- igraph
- MASS
- ggplot2

## Quick Start

```r
library(TRRident)

# Define markers
tcell_markers <- c('CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B2')
ex_markers <- c('PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'ENTPD1', 'CTLA4', 'TOX')

# 1. Load and prepare data
dat <- data_preparation_TRR(
  folder = 'path/to/visium/data',
  filenm = "filtered_feature_bc_matrix.h5",
  loc_file = 'path/to/spatial/tissue_positions.csv',
  res = 0.2,
  tcell_markers = tcell_markers,
  ex_markers = ex_markers
)

# 2. Estimate J1 parameter
j1_result <- estimate_J1_dat(genes = 'B2M', sp.dd = dat$sp.count.r)

# 3. Get lines and detect change points
all.lines <- get_lines(dat$loc.raw, grid.type = 'hex')

# 4. Run Gibbs sampling
ss.mat <- gibbs_sampling_new(...)

# 5. Get TRR regions
final.probs <- calculate_final_prob(ss.mat, ...)
final.regions <- get_TRR(final.probs)
```

## Method Overview

TRRident identifies T cell entry regions through the following steps:

1. **Data Preparation**: Load 10x Visium data and compute exhaustion scores adjusted for T cell marker expression
2. **J1 Estimation**: Estimate the coupling parameter J1 using variational inference on an Ising model
3. **Line Decomposition**: Decompose 2D spatial data into 1D lines along multiple directions
4. **Change Point Detection**: Detect change points in exhaustion marker expression using PELT with negative binomial regression
5. **J2 Estimation**: Estimate the impurity coupling parameter J2 at change points
6. **Gibbs Sampling**: Sample from the Ising model with impurities to get posterior probabilities
7. **Region Identification**: Identify connected regions with high TRR probability

## Citation

If you use TRRident in your research, please cite:

```
[Your publication details here]
```

## License

MIT License

## Contact

For questions or issues, please open an issue on GitHub or contact [your email].
