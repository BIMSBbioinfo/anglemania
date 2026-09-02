# =============================================================================
# Tests for Simpson's paradox management
# =============================================================================

# Helper: create a simple sparse matrix for testing
make_test_sparse <- function(n_genes = 50, n_cells = 200, seed = 42) {
    set.seed(seed)
    mat <- matrix(
        rpois(n_genes * n_cells, lambda = 2),
        nrow = n_genes, ncol = n_cells
    )
    rownames(mat) <- paste0("Gene", seq_len(n_genes))
    colnames(mat) <- paste0("Cell", seq_len(n_cells))
    Matrix::Matrix(mat, sparse = TRUE)
}

# =============================================================================
# 1. Microclustering tests
# =============================================================================

test_that("create_microclusters_sparse returns correct structure", {
    expr_mat <- make_test_sparse(n_genes = 50, n_cells = 200)
    result <- create_microclusters_sparse(
        expr_mat, n_clusters = 5, n_pcs = 5, min_cells = 10, verbose = FALSE
    )
    expect_type(result, "list")
    expect_named(result, c("assignments", "sizes", "n_clusters", "pca_coords"))
    expect_type(result$assignments, "integer")
    expect_true(result$n_clusters >= 1)
})

test_that("create_microclusters_sparse assigns all cells", {
    expr_mat <- make_test_sparse(n_genes = 50, n_cells = 200)
    result <- create_microclusters_sparse(
        expr_mat, n_clusters = 5, n_pcs = 5, min_cells = 10, verbose = FALSE
    )
    expect_equal(length(result$assignments), ncol(expr_mat))
    expect_true(all(names(result$assignments) == colnames(expr_mat)))
})

test_that("create_microclusters_sparse merges small clusters", {
    expr_mat <- make_test_sparse(n_genes = 50, n_cells = 200)
    result <- create_microclusters_sparse(
        expr_mat, n_clusters = 50, n_pcs = 5, min_cells = 30, verbose = FALSE
    )
    # With 200 cells and min_cells=30, can have at most 6 clusters
    expect_true(result$n_clusters <= 7)
    # All final clusters should have >= min_cells (30)
    expect_true(all(result$sizes >= 30))
})

test_that("create_microclusters_sparse limits clusters to cell count", {
    expr_mat <- make_test_sparse(n_genes = 50, n_cells = 60)
    result <- create_microclusters_sparse(
        expr_mat, n_clusters = 100, n_pcs = 5, min_cells = 20, verbose = FALSE
    )
    # With 60 cells and min_cells=20, at most 3 clusters
    expect_true(result$n_clusters <= 3)
})

test_that("create_microclusters_sparse uses irlba when available", {
    skip_if_not_installed("irlba")
    expr_mat <- make_test_sparse(n_genes = 100, n_cells = 300)
    result <- create_microclusters_sparse(
        expr_mat, n_clusters = 10, n_pcs = 10, min_cells = 20, verbose = FALSE
    )
    expect_true(result$n_clusters >= 1)
    expect_equal(sum(result$sizes), ncol(expr_mat))
})

test_that("create_microclusters_sparse warns when all clusters too small", {
    expr_mat <- make_test_sparse(n_genes = 50, n_cells = 20)
    # min_cells = 25 with only 20 cells — all clusters will be too small
    expect_warning(
        create_microclusters_sparse(
            expr_mat, n_clusters = 5, n_pcs = 3, min_cells = 25,
            verbose = FALSE
        ),
        "All clusters have fewer than"
    )
})

# =============================================================================
# 2. Weight computation tests
# =============================================================================

test_that("compute_simpsons_weights returns one weight per cell", {
    assignments <- c(A = 1L, B = 1L, C = 1L, D = 2L, E = 2L)
    w <- compute_simpsons_weights(assignments, min_cluster_size = 2)
    expect_equal(length(w), length(assignments))
    expect_true(all(names(w) == names(assignments)))
})

test_that("compute_simpsons_weights gives equal weight within cluster", {
    assignments <- c(
        A = 1L, B = 1L, C = 1L, D = 1L, E = 1L,
        F = 2L, G = 2L, H = 2L, I = 2L, J = 2L
    )
    w <- compute_simpsons_weights(assignments, min_cluster_size = 3)
    # All cells in cluster 1 should have same weight
    expect_equal(length(unique(w[assignments == 1])), 1)
    # All cells in cluster 2 should have same weight
    expect_equal(length(unique(w[assignments == 2])), 1)
    # Since both clusters have same size, weights should be equal
    expect_equal(w[["A"]], w[["F"]])
})

test_that("compute_simpsons_weights zeros small clusters", {
    assignments <- setNames(
        c(rep(1L, 20), rep(2L, 3)),
        paste0("Cell", seq_len(23))
    )
    w <- compute_simpsons_weights(assignments, min_cluster_size = 10)
    # Cluster 2 has only 3 cells → weight 0
    expect_true(all(w[assignments == 2] == 0))
    # Cluster 1 has 20 cells → weight > 0
    expect_true(all(w[assignments == 1] > 0))
})

test_that("compute_simpsons_weights falls back to uniform when all small", {
    assignments <- c(A = 1L, B = 1L, C = 2L, D = 2L)
    expect_warning(
        w <- compute_simpsons_weights(assignments, min_cluster_size = 10),
        "All clusters have fewer than"
    )
    expect_true(all(w > 0))
    expect_equal(length(unique(w)), 1)
})

test_that("compute_simpsons_weights: larger clusters get smaller per-cell weight", {
    assignments <- setNames(
        c(rep(1L, 30), rep(2L, 60)),
        paste0("Cell", seq_len(90))
    )
    w <- compute_simpsons_weights(assignments, min_cluster_size = 10)
    # Per-cell weight for cluster 1 should be larger than cluster 2
    expect_true(w[["Cell1"]] > w[["Cell31"]])
})

# =============================================================================
# 3. Weighted correlation tests
# =============================================================================

test_that("weighted_cor_fbm matches manual weighted Pearson for small matrix", {
    # 3 genes, 6 cells
    mat <- matrix(
        c(1, 2, 3, 4, 5, 6,
          2, 4, 6, 8, 10, 12,
          6, 5, 4, 3, 2, 1),
        nrow = 3, ncol = 6, byrow = TRUE
    )
    fbm <- bigstatsr::FBM(nrow = 3, ncol = 6, init = mat)
    w <- rep(1 / 6, 6)

    result <- weighted_cor_fbm(fbm, w)
    result_mat <- result[]

    # With uniform weights, should match standard Pearson
    expected <- cor(t(mat))
    diag(expected) <- NA

    expect_equal(result_mat, expected, tolerance = 1e-10)
})

test_that("weighted_cor_fbm with uniform weights matches extract_angles", {
    mat <- matrix(
        c(5, 3, 0, 0, 2,
          0, 0, 0, 3, 4,
          2, 1, 3, 4, 0,
          1, 2, 1, 2, 3),
        nrow = 4, ncol = 5, byrow = TRUE
    )
    fbm <- bigstatsr::FBM(nrow = 4, ncol = 5, init = mat)
    w <- rep(1, 5)

    result <- weighted_cor_fbm(fbm, w)[]

    # Standard correlation
    expected <- cor(t(mat))
    diag(expected) <- NA

    expect_equal(result, expected, tolerance = 1e-10)
})

test_that("weighted_cor_fbm is symmetric", {
    mat <- matrix(rpois(20, 3), nrow = 4, ncol = 5)
    fbm <- bigstatsr::FBM(nrow = 4, ncol = 5, init = mat)
    w <- c(1, 2, 1, 2, 1)

    result <- weighted_cor_fbm(fbm, w)[]
    # Check symmetry (ignoring NA diagonal)
    result_nodiag <- result
    diag(result_nodiag) <- 0
    expect_equal(result_nodiag, t(result_nodiag), tolerance = 1e-10)
})

test_that("weighted_cor_fbm has NA diagonal", {
    mat <- matrix(rpois(20, 3), nrow = 4, ncol = 5)
    fbm <- bigstatsr::FBM(nrow = 4, ncol = 5, init = mat)
    w <- rep(1, 5)

    result <- weighted_cor_fbm(fbm, w)[]
    expect_true(all(is.na(diag(result))))
})

test_that("weighted_cor_fbm output has correct dimensions", {
    n_genes <- 10
    n_cells <- 20
    mat <- matrix(rpois(n_genes * n_cells, 3), nrow = n_genes, ncol = n_cells)
    fbm <- bigstatsr::FBM(nrow = n_genes, ncol = n_cells, init = mat)
    w <- rep(1, n_cells)

    result <- weighted_cor_fbm(fbm, w)
    expect_equal(nrow(result), n_genes)
    expect_equal(ncol(result), n_genes)
})

test_that("weighted_cor_fbm handles zero-variance genes", {
    mat <- matrix(
        c(5, 5, 5, 5,
          1, 2, 3, 4,
          4, 3, 2, 1),
        nrow = 3, ncol = 4, byrow = TRUE
    )
    fbm <- bigstatsr::FBM(nrow = 3, ncol = 4, init = mat)
    w <- rep(1, 4)

    result <- weighted_cor_fbm(fbm, w)[]
    # Gene 1 has zero variance → correlations with it should be NA
    expect_true(all(is.na(result[1, ])))
    expect_true(all(is.na(result[, 1])))
})

test_that("weighted_cor_fbm values are in [-1, 1] or NA", {
    set.seed(99)
    mat <- matrix(rpois(50, 3), nrow = 5, ncol = 10)
    fbm <- bigstatsr::FBM(nrow = 5, ncol = 10, init = mat)
    w <- runif(10)

    result <- weighted_cor_fbm(fbm, w)[]
    valid <- result[!is.na(result)]
    expect_true(all(valid >= -1 - 1e-10 & valid <= 1 + 1e-10))
})

# =============================================================================
# 4. Local correlation tests
# =============================================================================

test_that("cor_from_stats returns NA for small n", {
    result <- cor_from_stats(10, 20, 30, 40, 50, n = 5, min_n = 20)
    expect_true(is.na(result))
})

test_that("cor_from_stats returns NA for zero-variance", {
    # All same values: var = 0
    result <- cor_from_stats(
        sum_x = 100, sum_y = 200, sum_xy = 2000,
        sum_x2 = 1000, sum_y2 = 4000, n = 10, min_n = 5
    )
    expect_true(is.na(result))
})

test_that("cor_from_stats matches manual computation", {
    x <- c(1, 2, 3, 4, 5)
    y <- c(2, 4, 6, 8, 10)
    n <- 5
    r <- cor_from_stats(
        sum_x = sum(x), sum_y = sum(y), sum_xy = sum(x * y),
        sum_x2 = sum(x^2), sum_y2 = sum(y^2), n = n, min_n = 3
    )
    expect_equal(r, cor(x, y), tolerance = 1e-10)
})

test_that("compute_local_correlations returns correct structure", {
    expr_mat <- make_test_sparse(n_genes = 20, n_cells = 100)
    assignments <- setNames(
        rep(1:2, each = 50),
        colnames(expr_mat)
    )
    gene_pairs <- data.frame(
        geneA = c("Gene1", "Gene2"),
        geneB = c("Gene3", "Gene4"),
        stringsAsFactors = FALSE
    )
    result <- compute_local_correlations(
        expr_mat, gene_pairs, assignments,
        min_cells_cor = 10, top_quantile = 0.95
    )
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), nrow(gene_pairs))
    expect_true(all(c("geneA", "geneB", "peak_cor", "peak_expr_level",
                       "n_valid_clusters") %in% colnames(result)))
})

test_that("compute_local_correlations handles missing genes", {
    expr_mat <- make_test_sparse(n_genes = 20, n_cells = 100)
    assignments <- setNames(
        rep(1:2, each = 50),
        colnames(expr_mat)
    )
    gene_pairs <- data.frame(
        geneA = c("Gene1", "NotAGene"),
        geneB = c("Gene2", "Gene3"),
        stringsAsFactors = FALSE
    )
    result <- compute_local_correlations(
        expr_mat, gene_pairs, assignments,
        min_cells_cor = 10, top_quantile = 0.95
    )
    expect_equal(nrow(result), 2)
    # "NotAGene" pair should have NA
    expect_true(is.na(result$peak_cor[2]))
})

# =============================================================================
# 5. Cross-dataset comparison tests
# =============================================================================

test_that("compare_local_correlations returns correct structure", {
    lc1 <- data.frame(
        geneA = c("G1", "G2"), geneB = c("G3", "G4"),
        peak_cor = c(0.8, 0.3), peak_expr_level = c(5, 2),
        n_valid_clusters = c(3, 2), stringsAsFactors = FALSE
    )
    lc2 <- data.frame(
        geneA = c("G1", "G2"), geneB = c("G3", "G4"),
        peak_cor = c(0.7, 0.1), peak_expr_level = c(4, 3),
        n_valid_clusters = c(3, 1), stringsAsFactors = FALSE
    )
    result <- compare_local_correlations(
        list(batch1 = lc1, batch2 = lc2),
        weights = c(batch1 = 1, batch2 = 1)
    )
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 2)
    expect_true(all(c("geneA", "geneB", "simpsons_score",
                       "mean_peak_cor", "sd_peak_cor") %in% colnames(result)))
})

test_that("compare_local_correlations: consistent pairs get higher score", {
    lc1 <- data.frame(
        geneA = c("G1", "G2"), geneB = c("G3", "G4"),
        peak_cor = c(0.8, 0.8), peak_expr_level = c(5, 5),
        n_valid_clusters = c(3, 3), stringsAsFactors = FALSE
    )
    lc2 <- data.frame(
        geneA = c("G1", "G2"), geneB = c("G3", "G4"),
        peak_cor = c(0.8, 0.1), peak_expr_level = c(5, 5),
        n_valid_clusters = c(3, 3), stringsAsFactors = FALSE
    )
    result <- compare_local_correlations(
        list(b1 = lc1, b2 = lc2),
        weights = c(b1 = 1, b2 = 1)
    )
    # G1-G3 is consistent (0.8, 0.8) → higher score
    # G2-G4 is inconsistent (0.8, 0.1) → lower score
    expect_true(result$simpsons_score[1] > result$simpsons_score[2])
})

test_that("compare_local_correlations handles single batch", {
    lc1 <- data.frame(
        geneA = "G1", geneB = "G2",
        peak_cor = 0.6, peak_expr_level = 3,
        n_valid_clusters = 2, stringsAsFactors = FALSE
    )
    result <- compare_local_correlations(
        list(batch1 = lc1),
        weights = c(batch1 = 1)
    )
    expect_equal(nrow(result), 1)
    expect_equal(result$simpsons_score, 0.6)
    expect_equal(result$sd_peak_cor, 0)
})

test_that("compare_local_correlations handles NA peak_cor", {
    lc1 <- data.frame(
        geneA = "G1", geneB = "G2",
        peak_cor = NA_real_, peak_expr_level = NA_real_,
        n_valid_clusters = 0, stringsAsFactors = FALSE
    )
    lc2 <- data.frame(
        geneA = "G1", geneB = "G2",
        peak_cor = 0.5, peak_expr_level = 3,
        n_valid_clusters = 2, stringsAsFactors = FALSE
    )
    result <- compare_local_correlations(
        list(b1 = lc1, b2 = lc2),
        weights = c(b1 = 1, b2 = 1)
    )
    # Only batch2 has valid data → score = mean_peak_cor from single batch
    expect_equal(result$simpsons_score, 0.5)
})

# =============================================================================
# 6. Pipeline orchestration tests
# =============================================================================

test_that("prepare_simpsons returns correct structure", {
    mat1 <- make_test_sparse(n_genes = 50, n_cells = 200, seed = 1)
    mat2 <- make_test_sparse(n_genes = 50, n_cells = 150, seed = 2)
    matrix_list <- list(batch1 = mat1, batch2 = mat2)

    result <- prepare_simpsons(
        matrix_list,
        n_clusters = 5, n_pcs = 5, min_cells = 20,
        min_cluster_size = 15, verbose = FALSE
    )
    expect_type(result, "list")
    expect_named(result, c("weights_per_batch", "clusters_per_batch"))
    expect_named(result$weights_per_batch, c("batch1", "batch2"))
    expect_named(result$clusters_per_batch, c("batch1", "batch2"))
    # Weights have correct length
    expect_equal(length(result$weights_per_batch$batch1), 200)
    expect_equal(length(result$weights_per_batch$batch2), 150)
})

test_that("factorise with cell_weights runs without error", {
    mat <- matrix(rpois(30, 3), nrow = 6, ncol = 5)
    fbm <- bigstatsr::FBM(nrow = 6, ncol = 5, init = mat)
    w <- rep(1, 5)

    result <- factorise(fbm, method = "cosine", seed = 1, cell_weights = w)
    expect_s4_class(result, "FBM")
    expect_equal(nrow(result), 6)
    expect_equal(ncol(result), 6)
})

test_that("factorise with uniform cell_weights gives same shape as standard", {
    mat <- matrix(
        c(5, 3, 0, 0, 2,
          0, 0, 0, 3, 4,
          2, 1, 3, 4, 0,
          1, 2, 1, 2, 3),
        nrow = 4, ncol = 5, byrow = TRUE
    )
    fbm1 <- bigstatsr::FBM(nrow = 4, ncol = 5, init = mat)
    fbm2 <- bigstatsr::FBM(nrow = 4, ncol = 5, init = mat)

    result_std <- factorise(fbm1, method = "cosine", seed = 1)
    result_w <- factorise(
        fbm2, method = "cosine", seed = 1,
        cell_weights = rep(1, 5)
    )
    # Same dimensions
    expect_equal(dim(result_std), dim(result_w))
    # Both should be gene x gene z-scores
    expect_equal(nrow(result_std), 4)
    expect_equal(ncol(result_std), 4)
})

test_that("factorise with spearman + cell_weights warns", {
    mat <- matrix(rpois(30, 3), nrow = 6, ncol = 5)
    fbm <- bigstatsr::FBM(nrow = 6, ncol = 5, init = mat)
    w <- rep(1, 5)

    expect_warning(
        factorise(fbm, method = "spearman", seed = 1, cell_weights = w),
        "Spearman not supported with Simpson's correction"
    )
})

# =============================================================================
# 7. check_params integration tests
# =============================================================================

test_that("check_params accepts Simpson's parameters", {
    sce <- sce_example()
    params <- check_params(
        sce,
        batch_key = "batch",
        dataset_key = "dataset",
        max_n_genes = 2000,
        method = "cosine",
        min_cells_per_gene = 1,
        min_samples_per_gene = 2,
        allow_missing_features = FALSE,
        permute_row_or_column = "column",
        permutation_function = "sample",
        prefilter_threshold = 0.5,
        normalization_method = "divide_by_total_counts",
        verbose = FALSE,
        use_simpsons = TRUE,
        simpsons_n_clusters = 50,
        simpsons_n_pcs = 10
    )
    expect_true(params$use_simpsons)
    expect_equal(params$simpsons_n_clusters, 50)
    expect_equal(params$simpsons_n_pcs, 10)
    # Defaults
    expect_equal(params$simpsons_min_cells, 30)
    expect_equal(params$simpsons_weight_in_ranking, 0.0)
})

test_that("check_params defaults don't break without Simpson's args", {
    sce <- sce_example()
    params <- check_params(
        sce,
        batch_key = "batch",
        dataset_key = "dataset",
        max_n_genes = 2000,
        method = "cosine",
        min_cells_per_gene = 1,
        min_samples_per_gene = 2,
        allow_missing_features = FALSE,
        permute_row_or_column = "column",
        permutation_function = "sample",
        prefilter_threshold = 0.5,
        normalization_method = "divide_by_total_counts",
        verbose = FALSE
    )
    expect_false(params$use_simpsons)
})

test_that("check_params validates Simpson's params when enabled", {
    sce <- sce_example()
    expect_error(
        check_params(
            sce,
            batch_key = "batch",
            dataset_key = "dataset",
            max_n_genes = 2000,
            method = "cosine",
            min_cells_per_gene = 1,
            min_samples_per_gene = 2,
            allow_missing_features = FALSE,
            permute_row_or_column = "column",
            permutation_function = "sample",
            prefilter_threshold = 0.5,
            normalization_method = "divide_by_total_counts",
            verbose = FALSE,
            use_simpsons = TRUE,
            simpsons_n_clusters = -1
        )
    )
})

# =============================================================================
# 8. Full pipeline integration test
# =============================================================================

test_that("anglemania with use_simpsons=TRUE runs on sce_example", {
    set.seed(1)
    sce <- sce_example()
    # Basic Simpson's correction (reweighting only, no posthoc)
    suppressMessages(
        sce_result <- anglemania(
            sce,
            batch_key = "batch",
            method = "cosine",
            use_simpsons = TRUE,
            simpsons_n_clusters = 3,
            simpsons_n_pcs = 5,
            simpsons_min_cells = 10,
            simpsons_min_cluster_size = 5,
            verbose = FALSE
        )
    )
    # Check that simpsons_data is in metadata
    expect_true(
        !is.null(S4Vectors::metadata(sce_result)$anglemania$simpsons_data)
    )
    # Check that genes were selected
    genes <- get_anglemania_genes(sce_result)
    expect_true(length(genes) > 0)
})

test_that("anglemania with use_simpsons=FALSE is unchanged", {
    set.seed(1)
    sce <- sce_example()
    suppressMessages(
        sce_result <- anglemania(
            sce,
            batch_key = "batch",
            method = "cosine",
            use_simpsons = FALSE,
            verbose = FALSE
        )
    )
    # No simpsons_data
    expect_null(
        S4Vectors::metadata(sce_result)$anglemania$simpsons_data
    )
    genes <- get_anglemania_genes(sce_result)
    expect_true(length(genes) > 0)
})

# =============================================================================
# Binary PCA tests
# =============================================================================

test_that("create_microclusters_sparse with binary PCA returns valid result", {
    expr_mat <- make_test_sparse(n_genes = 50, n_cells = 200)
    result <- create_microclusters_sparse(
        expr_mat, n_clusters = 5, n_pcs = 5, min_cells = 10,
        use_binary_pca = TRUE, verbose = FALSE
    )
    expect_type(result, "list")
    expect_named(result, c("assignments", "sizes", "n_clusters", "pca_coords"))
    expect_type(result$assignments, "integer")
    expect_true(result$n_clusters >= 1)
    expect_equal(nrow(result$pca_coords), 200)
})

test_that("binary PCA produces different clusters than expression PCA", {
    expr_mat <- make_test_sparse(n_genes = 50, n_cells = 200)
    res_expr <- create_microclusters_sparse(
        expr_mat, n_clusters = 5, n_pcs = 5, min_cells = 10,
        use_binary_pca = FALSE, verbose = FALSE
    )
    res_bin <- create_microclusters_sparse(
        expr_mat, n_clusters = 5, n_pcs = 5, min_cells = 10,
        use_binary_pca = TRUE, verbose = FALSE
    )
    # PCA coordinates should differ
    expect_false(identical(res_expr$pca_coords, res_bin$pca_coords))
})

# =============================================================================
# Density weights tests
# =============================================================================

test_that("compute_density_weights returns valid weights", {
    set.seed(42)
    pca <- matrix(rnorm(200), nrow = 100, ncol = 2)
    rownames(pca) <- paste0("Cell", 1:100)
    w <- compute_density_weights(pca, k = 10)
    expect_length(w, 100)
    expect_true(all(w > 0))
    expect_true(all(is.finite(w)))
    # Weights should sum to n_cells
    expect_equal(sum(w), 100, tolerance = 0.01)
})

test_that("density weights upweight sparse regions", {
    set.seed(42)
    # Create two clusters: one dense, one sparse
    dense <- matrix(rnorm(160, mean = 0, sd = 0.1), nrow = 80, ncol = 2)
    sparse <- matrix(rnorm(40, mean = 5, sd = 2), nrow = 20, ncol = 2)
    pca <- rbind(dense, sparse)
    rownames(pca) <- paste0("Cell", 1:100)
    w <- compute_density_weights(pca, k = 10)
    # Sparse region cells should have higher weights
    expect_true(mean(w[81:100]) > mean(w[1:80]))
})

test_that("prepare_simpsons with density weights runs", {
    expr_mat <- make_test_sparse(n_genes = 50, n_cells = 200)
    result <- prepare_simpsons(
        matrix_list = list(batch1 = expr_mat),
        n_clusters = 5, n_pcs = 5, min_cells = 10,
        use_density_weights = TRUE, density_k = 10,
        verbose = FALSE
    )
    expect_true(length(result$weights_per_batch$batch1) == 200)
    expect_true(all(result$weights_per_batch$batch1 > 0))
})

test_that("prepare_simpsons with binary PCA + density weights runs", {
    expr_mat <- make_test_sparse(n_genes = 50, n_cells = 200)
    result <- prepare_simpsons(
        matrix_list = list(batch1 = expr_mat),
        n_clusters = 5, n_pcs = 5, min_cells = 10,
        use_binary_pca = TRUE,
        use_density_weights = TRUE, density_k = 10,
        verbose = FALSE
    )
    expect_true(length(result$weights_per_batch$batch1) == 200)
    expect_true(all(result$weights_per_batch$batch1 > 0))
})
