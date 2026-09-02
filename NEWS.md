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

* `simpsons_n_clusters` now defaults to `NULL`, meaning the natural community
structure Leiden finds at resolution 1 — which is what every published benchmark
actually used, and the best-performing setting. Previously it defaulted to `15`
but the resolution was computed as `max(1, simpsons_n_clusters / 20)` and floored
at 1, so every value from 2 to 20 gave identical clustering and the argument had
no effect. Supplying an integer now genuinely targets that many clusters, by
bisecting on resolution (the resolution needed depends on the dataset and its
size, so no fixed formula can work).

* Forcing finer clusterings is worse, not better. On SCIB Sim2, DE precision is
0.983 with `NULL`, against 0.903 for no correction at all, 0.812 at 15 clusters,
0.800 at 30 and 0.857 at 60. The correction works by giving each cell type equal
weight; splitting cell types into finer microclusters removes the between-type
variation that carries the correlation signal. Despite the name, this is
cell-type reweighting rather than microclustering.

* The correction needs enough cells per cluster: cluster-inverse weights promote
small clusters to full weight while they remain noisy, so batches with small or
highly skewed populations can come out worse than uncorrected.

# anglemania 0.99.5

* Initial submission to Bioconductor.

* anglemania provides improved feature extraction for scRNA-seq dataset integration.

* Checkout the vignette for a step-by-step tutorial on how to use anglemania. 

* anglemania can be used on top of `r Biocpkg("SingleCellExperiment")` or 
`r Biocpkg("SummarizedExperiment")` objects.