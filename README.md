# Bayesian Ordinal Factor MFM Clustering

Code for **“Bayesian Ordinal Factor Mixture of Finite Mixtures Modeling for
Single-Cell DNA Methylation Clustering.”**

## Folder structure

```text
function/
  ordinal_mfm.R     Functions called by the user
  mcmc.cpp          MCMC implementation

data/
  simulated_balanced_strong_n1000_p5000.rds

example/
  run_example.R     Complete example analysis
```

## Example data

The included dataset follows the model-based **balanced, strong-signal**
simulation setting used in the paper:

- 5,000 methylation features;
- 1,000 cells;
- five balanced true clusters of 200 cells each;
- strong-signal parameters `s = 5`, `p_sig = 0.40`, and `sigma_xi = 0.20`;
- anisotropy `8` and cluster volume multipliers
  `(0.6, 0.8, 1.0, 1.4, 2.0)`;
- generation seed `20009`.

Load it in R:

```r
example_data <- readRDS("data/simulated_balanced_strong_n1000_p5000.rds")
Y <- example_data$Y                    # 5000 features × 1000 cells
true_label <- example_data$true_label  # used only for evaluation
example_data$metadata                  # simulation settings
```

The file also contains `W`, the corresponding continuous methylation
proportion matrix. The model is fitted to the five-category ordinal matrix `Y`.

## Install required packages

```r
install.packages(c("Rcpp", "RcppArmadillo", "mclust"))
```

A working C++ compiler for R is required.

## Use the function

Run R from the folder containing this README:

```r
source("function/ordinal_mfm.R")

fit <- fit_ordinal_mfm(Y, seed = 123)
result <- summarize_clustering(fit)
```

The default MCMC configuration is:

| Setting | Default |
|---|---:|
| Latent factor dimension `q` | 10 |
| Independent MCMC chains | 4 |
| Iterations per chain | 2,000 |
| Burn-in per chain | 1,000 |
| Initial clusters | 10 |
| `sigma_B` | 5 |
| `alpha_mfm` | 1 |
| Maximum occupied clusters | 40 |
| `kappa0` | 0.1 |
| `nu0` | `q + 2` |
| `lambda_pois` | 9 |

The four chains use consecutive seeds beginning with `seed`. For example,
`seed = 123` uses seeds 123, 124, 125, and 126. Post-burn-in partitions from
all chains are combined for the posterior similarity matrix and Dahl partition.
The original results for every chain remain available in `fit$chains`.

Main output:

```r
result$cluster                    # Dahl label for every cell
result$n_clusters                 # clusters in the Dahl partition
result$posterior_mode_n_clusters
result$posterior_mode_by_chain    # separate diagnostic for each chain
result$posterior_similarity       # 1000 × 1000 co-clustering matrix
result$ncluster_trace_by_chain
```

`q` is the factor dimension, not a prespecified cluster count. The MFM model
infers the number of occupied clusters.

## Run the complete example

From Terminal:

```bash
cd /path/to/this/folder
Rscript example/run_example.R
```

This runs all four full MCMC chains. With 5,000 features and 1,000 cells, it is
a substantive analysis rather than a quick smoke test and may take considerable
time and memory.

When the run finishes, the following files are saved in `example/`:

- `example_metrics.csv`;
- `example_cluster_assignments.csv`;
- `example_results.rds`, containing all four chains and the clustering summary.

The output files are created only after a successful run.

## Input for a new analysis

`Y` must have features in rows, cells in columns, values in `{1,2,3,4,5}`, and
no missing values. Preprocessing and imputation must be completed first.

For a complete methylation proportion matrix `W` in `[0,1]`:

```r
Y <- methylation_to_ordinal(W)
fit <- fit_ordinal_mfm(Y, seed = 123)
result <- summarize_clustering(fit)
```

The five intervals are `[0,0.2]`, `(0.2,0.4]`, `(0.4,0.6]`, `(0.6,0.8]`,
and `(0.8,1]`.
