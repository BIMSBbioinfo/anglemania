# anglemania 0.99.7

* Added Simpson's paradox correction, enabled with `anglemania(use_simpsons = TRUE)`.
Cells are grouped into microclusters within each batch (PCA + SNN graph + Leiden) and
gene-gene correlations are computed with per-cell weights `1 / |cluster|`, so every
microcluster contributes equally regardless of size. This removes correlation signal
driven purely by differing cell type proportions across batches. On SCIB simulations
DE precision rises to 0.937-0.983 from a baseline of 0.899-0.901.

* Added `simpsons_use_count_split`, which uses count splitting
([Neufeld et al.](https://arxiv.org/abs/2307.12985)) so that clustering and correlation
estimation run on independent Poisson folds of the counts. Recommended on real data,
where uncorrected Simpson's weighting tends to overcorrect. Requires the `countsplit`
package.

* Added `simpsons_weight_in_ranking`, an optional third criterion in `select_genes()`
that scores gene pairs on how consistent their within-microcluster correlations are
across datasets.

* Fixed `permute_nonzero()` erroring on vectors with fewer than two non-zero entries.

# anglemania 0.99.5

* Initial submission to Bioconductor.

* anglemania provides improved feature extraction for scRNA-seq dataset integration.

* Checkout the vignette for a step-by-step tutorial on how to use anglemania. 

* anglemania can be used on top of `r Biocpkg("SingleCellExperiment")` or 
`r Biocpkg("SummarizedExperiment")` objects.