# anglemania 0.99.7

* Added Simpson's paradox correction, enabled with `anglemania(use_simpsons = TRUE)`.
Cells are grouped into clusters within each batch (PCA + SNN graph + Leiden) and
gene-gene correlations are computed with per-cell weights `1 / |cluster|`, so every
cluster contributes equally regardless of size. This removes correlation signal
driven purely by differing cell type proportions across batches. On SCIB simulations
DE precision rises to 0.937-0.983 from a baseline of 0.899-0.901.

* Added `simpsons_use_count_split`, which uses count splitting
([Neufeld et al.](https://arxiv.org/abs/2307.12985)) so that clustering and correlation
estimation run on independent Poisson folds of the counts. Recommended on real data,
where uncorrected Simpson's weighting tends to overcorrect. Requires the `countsplit`
package.

* Added `simpsons_weight_in_ranking`, an optional third criterion in `select_genes()`
that scores gene pairs on how consistent their within-cluster correlations are
across datasets.

* Exported `create_microclusters_sparse()` and `compute_simpsons_weights()` so the
clustering and the per-cell weights can be inspected directly.

* Known limitation: `simpsons_n_clusters` does not control clustering granularity.
It sets the Leiden resolution to `max(1, simpsons_n_clusters / 20)`, so every value
from 2 to 20 gives resolution 1 and identical clustering, and above 20
`igraph::cluster_leiden(objective_function = "modularity")` responds only weakly.
The clustering resolves whole cell types rather than finer sub-structures. The
benchmark results above were produced at the default and are unaffected, but the
argument cannot currently be tuned. The correction also needs enough cells per
cluster: cluster-inverse weights promote small clusters to full weight while they
remain noisy, so batches with small or highly skewed populations can come out worse.

* Fixed `permute_nonzero()` erroring on vectors with fewer than two non-zero entries.

# anglemania 0.99.5

* Initial submission to Bioconductor.

* anglemania provides improved feature extraction for scRNA-seq dataset integration.

* Checkout the vignette for a step-by-step tutorial on how to use anglemania. 

* anglemania can be used on top of `r Biocpkg("SingleCellExperiment")` or 
`r Biocpkg("SummarizedExperiment")` objects.