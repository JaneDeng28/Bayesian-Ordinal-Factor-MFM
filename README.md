# SMORE: Single-cell MethylOme Reduction and Embedding

This repository provides the R and C++ implementation of **SMORE** (**S**ingle-cell **M**ethyl**O**me **R**eduction and **E**mbedding), a Bayesian framework for identifying cell populations from single-cell DNA methylation data. SMORE analyzes a preprocessed ordinal methylation matrix, with genomic regions represented by rows and individual cells by columns, and jointly learns a low-dimensional representation, cell-population assignments, and the number of occupied populations.

SMORE links observed methylation states to latent continuous variables through an ordinal probit model. A low-rank factor model captures coordinated methylation variation across genomic regions, while a mixture-of-finite-mixtures (MFM) prior clusters cells through their latent factor scores. Because the observation model, latent representation, and clustering structure are estimated jointly, uncertainty can propagate throughout the analysis rather than being separated across independent dimension-reduction and clustering steps.

Posterior samples are summarized using a cell-by-cell co-clustering probability matrix and Dahl's least-squares representative partition. The resulting output includes representative cell labels, posterior inference on the number of occupied populations, pairwise clustering uncertainty, and chain-specific diagnostic summaries.

## Workflow

1.  Prepare a feature-by-cell ordinal methylation matrix with values 1–5, or convert a complete methylation proportion matrix with `methylation_to_ordinal()`.
2.  Fit the model with `fit_ordinal_mfm()`.
3.  Summarize the posterior clustering with `summarize_clustering()`.
4.  Inspect the representative labels, posterior cluster count, co-clustering probabilities, and chain-specific diagnostics.

## Installation

Install the required R packages:

``` r
install.packages(c("Rcpp", "RcppArmadillo", "mclust"))
```

A working C++ compiler for R is required. `Rcpp` and `RcppArmadillo` are needed to fit the model; `mclust` is used only to calculate the adjusted Rand index in the simulated-data example.

Clone or download this repository, start R in the repository root, and load the user-facing functions:

``` r
source("function/ordinal_mfm.R")
```

The C++ sampler is compiled automatically on the first call to `fit_ordinal_mfm()`.

## Input data

### Ordinal input

The main function accepts `Y`, a numeric or integer matrix with:

| Property       | Required format                               |
|----------------|-----------------------------------------------|
| Rows           | Methylation features or genomic regions       |
| Columns        | Cells                                         |
| Entries        | Ordinal categories `1`, `2`, `3`, `4`, or `5` |
| Missing values | Not allowed                                   |
| Minimum size   | At least one feature and two cells            |

For example, an input with four features and three cells looks like this:

``` r
Y <- matrix(
  c(1, 2, 1,
    4, 5, 4,
    2, 3, 2,
    5, 4, 5),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(
    paste0("Feature_", 1:4),
    paste0("Cell_", 1:3)
  )
)

Y
#           Cell_1 Cell_2 Cell_3
# Feature_1      1      2      1
# Feature_2      4      5      4
# Feature_3      2      3      2
# Feature_4      5      4      5
```

Cell names are optional. The order of `result$cluster` follows the column order of `Y`.

The current implementation is written specifically for five ordinal levels. This is reflected in the R input validation, the five-category cutpoint vector, and the C++ sampler. Therefore, the number of ordinal levels is not currently a user-configurable argument: changing from five levels to another number would require corresponding code changes, not merely a different value of `q`.

### Continuous methylation proportions

If the starting data are a complete feature-by-cell methylation proportion matrix `W` in `[0, 1]`, convert it before fitting:

``` r
Y <- methylation_to_ordinal(W)
```

The conversion uses the intervals `[0, 0.2]`, `(0.2, 0.4]`, `(0.4, 0.6]`, `(0.6, 0.8]`, and `(0.8, 1]`, mapped to categories 1–5. Equivalently, use the convenience wrapper:

``` r
fit <- fit_methylation_mfm(W, seed = 123)
```

Both `Y` and `W` must be complete: preprocessing and imputation must be performed before calling these functions.

## Quick start

``` r
source("function/ordinal_mfm.R")

example_data <- readRDS("data/example_data_n1000_p5000.rds")
Y <- example_data$Y

fit <- fit_ordinal_mfm(Y, q = 4, seed = 123)
result <- summarize_clustering(fit)

head(result$cluster)
result$n_clusters
result$posterior_mode_n_clusters
result$posterior_mode_by_chain
```

`true_label` is included in the simulated example for evaluation only. It is not supplied to the model.

## Main functions

### `fit_ordinal_mfm()`

``` r
fit <- fit_ordinal_mfm(
  Y,
  q = 4,
  chains = 4,
  n_iter = 2000,
  burn_in = 1000,
  initial_clusters = 10,
  seed = 123
)
```

Important arguments are:

| Argument | Description | Default |
|----|----|---:|
| `Y` | Feature-by-cell ordinal matrix with entries in `{1, ..., 5}` | required |
| `q` | Dimension of each cell's latent factor score; this is neither the number of ordinal levels nor the number of clusters | `4` |
| `chains` | Number of independent MCMC chains | `4` |
| `n_iter` | Total iterations per chain | `2000` |
| `burn_in` | Initial iterations discarded from each chain | `1000` |
| `initial_clusters` | Cluster count used only to initialize the sampler | `10` |
| `sigma_B` | Prior standard deviation for free factor loadings | `5` |
| `alpha_mfm` | MFM Dirichlet concentration parameter | `1` |
| `max_clusters` | Computational upper bound on occupied clusters | `40` |
| `kappa0` | Normal-inverse-Wishart mean-precision parameter | `0.1` |
| `nu0_add` | Sets inverse-Wishart degrees of freedom to `q + nu0_add` | `2` |
| `lambda_pois` | Poisson prior mean for the component count minus one | `9` |
| `seed` | Base random seed; chain `c` uses `seed + c - 1` | `NULL` |

The return value, `fit`, is an object of class `ordinal_mfm_fit`:

| Component | Contents |
|----|----|
| `fit$chains` | Raw posterior output for every MCMC chain |
| `fit$input` | Number of features, number of cells, and factor dimension `q` |
| `fit$settings` | MCMC settings, hyperparameters, and chain seeds |

Chains are run sequentially to avoid multiplying peak memory use. When `seed = 123` and `chains = 4`, the chain seeds are 123, 124, 125, and 126.

#### What `q` means in the R and C++ code

In the manuscript, `q` is the dimension of the latent factor score `eta_j`. Accordingly, the loading matrix has dimensions `p × q`, each cell has a `q`-dimensional factor score, and the MFM clusters cells in this `q`-dimensional space. The current R interface defaults to `q = 4`, matching the latent dimension used to generate the included model-based example.

The C++ sampler names this same argument `K`, which can be confused with the manuscript notation for the number of mixture components. In the implementation, however, the C++ argument `K` is used to set the number of columns of the factor loading matrix `B` and the number of rows of the factor-score matrix `eta`:

``` text
R interface:       q
C++ sampler:       K        (the same latent factor dimension)
Loading matrix:    B        p × q
Factor scores:     eta      q × n
Initial clusters:  K_init   initialization only
Occupied clusters: n_clust  inferred by the MFM sampler
```

Thus, changing `q` changes the dimension of the latent representation; it does not prescribe how many cell clusters the model must return. The argument `initial_clusters` is passed to C++ as `K_init` and affects only the starting partition. The posterior number of occupied clusters is learned by the MFM model and is summarized by `n_clusters` and `posterior_mode_n_clusters`.

The latent dimension `q` is also separate from the number of ordinal categories in `Y`. Five ordinal categories describe the measurement scale of each matrix entry, whereas `q` describes the dimension of each cell's latent coordinate. There is no general rule that a model with five ordinal levels must use `q = 4`.

### `summarize_clustering()`

``` r
result <- summarize_clustering(fit)
```

This function combines all retained allocation draws, calculates the posterior similarity matrix, and selects Dahl's least-squares representative partition. It returns:

| Component | Shape | Meaning |
|----|---:|----|
| `cluster` | `n_cells` | Dahl cluster label for each input cell |
| `n_clusters` | 1 | Number of clusters in the Dahl partition |
| `posterior_mode_n_clusters` | 1 | Most frequent occupied-cluster count across all retained draws |
| `posterior_similarity` | `n_cells × n_cells` | Posterior probability that each pair of cells co-clusters |
| `ncluster_trace` | total retained draws | Occupied-cluster count across all chains |
| `ncluster_trace_by_chain` | list of length `chains` | Cluster-count trace for each chain |
| `posterior_mode_by_chain` | `chains` | Modal cluster count for each chain |
| `dahl_loss` | total retained draws | Dahl loss for every candidate partition |
| `selected_draw` | 1 | Index of the allocation draw selected by Dahl's criterion |

The numeric cluster labels are identifiers only: label 1 is not intrinsically greater than label 2. For scientific interpretation, investigate the methylation profiles of the cells assigned to each cluster. Large differences in `posterior_mode_by_chain` or unstable cluster-count traces should be examined before interpreting the representative partition.

## Included example data and output

`data/example_data_n1000_p5000.rds` follows the model-based balanced, strong-signal simulation setting used in the paper. It contains:

| Object       | Description                                                  |
|--------------|--------------------------------------------------------------|
| `Y`          | `5000 × 1000` ordinal input matrix with values 1–5           |
| `W`          | Corresponding continuous methylation proportions in `[0, 1]` |
| `true_label` | Five balanced simulated clusters of 200 cells each           |
| `metadata`   | Complete data-generation settings                            |

The simulation uses a four-dimensional latent signal (`q_signal = 4`), `s = 5`, `p_signal = 0.40`, `sigma_xi = 0.20`, anisotropy 8, volume multipliers `(0.6, 0.8, 1.0, 1.4, 2.0)`, and generation seed 20009.

### Why the example uses four latent dimensions

The values `q_signal = 4` and five ordinal categories arise from different parts of the simulation and should not be conflated. The example contains five true cell clusters. Their latent centers were constructed as the five vertices of a regular simplex so that all pairs of cluster centers were equally distant. A regular simplex with five vertices has affine dimension `5 - 1 = 4`. Therefore, four is the smallest latent dimension in which five equidistant cluster centers can be represented without distortion.

This geometric construction explains the data-generating choice `q_signal = K_true - 1 = 4` for `K_true = 5`. It was used to create five symmetrically separated clusters with the minimum required latent dimension. It does **not** imply that fitting with `q = 4` fixes the estimated number of clusters at five: the MFM sampler still infers the occupied cluster count.

Separately, the continuous methylation values were discretized into five ordinal levels. The manuscript uses five levels as a compromise between a binary methylated/unmethylated representation and overly fine discretization of noisy methylation proportions. This five-level measurement choice is not the mathematical reason for the four-dimensional simplex.

The checked-in example output was generated previously with `q = 10` and inferred five clusters:

``` text
metric                         value
dahl_clusters                 5
posterior_mode_clusters       5
adjusted_rand_index           0.943834450256832

posterior_mode_by_chain
chain_1  chain_2  chain_3  chain_4
      5        5        5        5
```

The adjusted Rand index compares the estimated labels with the known simulated labels. It is available for this benchmark because `true_label` is known; it cannot generally be calculated for an unlabeled real dataset.

These saved results document that earlier run and are not the output of the current `q = 4` default. Re-running `example/run_example.R` with the current R interface uses `q = 4` because the script does not explicitly override `q`.

## Reproduce the complete example

From a terminal in the repository root:

``` bash
Rscript example/run_example.R
```

With the current interface, this runs four full MCMC chains with `q = 4`, 2,000 iterations per chain, 5,000 features, and 1,000 cells. It is a substantive analysis and may require considerable time and memory. Files are written only after a successful run:

- `example/example_metrics.csv`: inferred cluster counts and adjusted Rand index;
- `example/example_cluster_assignments.csv`: cell IDs, simulated truth, and estimated labels;
- `example/example_results.rds`: all four raw chains and the combined clustering summary.

## Repository structure

``` text
function/
  ordinal_mfm.R       User-facing R functions
  mcmc.cpp            C++ MCMC implementation

data/
  example_data_n1000_p5000.rds

example/
  run_example.R
  example_metrics.csv
  example_cluster_assignments.csv
  example_results.rds
```
