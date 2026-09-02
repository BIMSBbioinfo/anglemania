# =============================================================================
# Simpson's Paradox Management for Anglemania
#
# Corrects composition-driven spurious correlations by:
#   1. Cluster-reweighting: each microcluster contributes equally to
#      correlation estimates, regardless of cluster size.
#   2. Local correlation preservation: computes per-microcluster correlations
#      and compares peaks across datasets to identify truly preserved circuits.
# =============================================================================


# =============================================================================
# PART A: Microclustering and cluster-reweighted correlations
# =============================================================================

#' Create microclusters from a sparse expression matrix
#'
#' Performs PCA on log-normalized expression, builds an SNN graph on the
#' PCA embedding and clusters cells into microclusters with Leiden. Small
#' clusters are merged into their nearest larger neighbor.
#'
#' @param expr_mat Sparse matrix (genes x cells), raw counts.
#' @param n_clusters Nominal target number of microclusters. Sets the Leiden
#'   resolution to \code{max(1, n_clusters / 20)}; because of the floor, all
#'   values <= 20 give resolution 1 and identical clustering. Default
#'   \code{15}.
#' @param n_pcs Number of principal components. Default \code{20}.
#' @param min_cells Minimum cells per microcluster; smaller clusters
#'   are merged. Default \code{30}.
#' @param use_binary_pca Logical; build the embedding from the binarised
#'   detection matrix instead of normalized expression. Default \code{FALSE}.
#' @param verbose Logical; print progress? Default \code{TRUE}.
#' @return A list with components:
#'   \describe{
#'     \item{assignments}{Named integer vector of cluster assignments.}
#'     \item{sizes}{Table of cluster sizes.}
#'     \item{n_clusters}{Final number of clusters.}
#'   }
#'
#' @importFrom Matrix colSums rowMeans
#' @importFrom stats prcomp median var
#'
#' @examples
#' sce <- sce_example()
#' mat <- SingleCellExperiment::counts(sce)
#' clusters <- create_microclusters_sparse(mat, n_clusters = 15, min_cells = 20)
#' clusters$n_clusters
#' table(clusters$assignments)
#'
#' @export
create_microclusters_sparse <- function(
    expr_mat,
    n_clusters = 15,
    n_pcs = 20,
    min_cells = 30,
    use_binary_pca = FALSE,
    verbose = TRUE
) {
    required <- c("bluster", "igraph")
    missing_pkgs <- required[
        !vapply(required, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_pkgs) > 0) {
        stop(
            "Simpson's correction requires the package(s): ",
            paste(missing_pkgs, collapse = ", "),
            ". Install with: BiocManager::install(c(",
            paste0("\"", missing_pkgs, "\"", collapse = ", "),
            "))"
        )
    }

    n_cells <- ncol(expr_mat)

    if (verbose) {
        message(sprintf(
            "  Creating %s microclusters from %s cells%s...",
            n_clusters, n_cells,
            if (use_binary_pca) " (binary PCA)" else ""
        ))
    }

    if (use_binary_pca) {
        # Binary PCA: use detection matrix (0/1) instead of expression
        # This avoids using gene expression values for manifold construction,
        # reducing circularity with downstream correlation computation
        binary_mat <- (expr_mat > 0) + 0

        # Subsample genes for speed
        if (nrow(binary_mat) > 2000) {
            # Select genes with intermediate detection rates (most informative)
            det_rate <- Matrix::rowMeans(binary_mat)
            informative <- det_rate > 0.01 & det_rate < 0.99
            if (sum(informative) > 2000) {
                gv <- det_rate[informative] * (1 - det_rate[informative])
                top_genes <- names(sort(gv, decreasing = TRUE))[seq_len(2000)]
            } else {
                top_genes <- names(which(informative))
            }
            binary_sub <- binary_mat[top_genes, ]
        } else {
            binary_sub <- binary_mat
        }

        # Centre by gene detection rates
        gene_means <- Matrix::rowMeans(binary_sub)
        centered <- binary_sub - gene_means

        # PCA via truncated SVD
        n_pcs_use <- min(n_pcs,
                         min(nrow(centered), ncol(centered)) - 1)
        if (requireNamespace("irlba", quietly = TRUE)) {
            svd_res <- irlba::irlba(Matrix::t(centered), nv = n_pcs_use)
            pca_coords <- svd_res$u[, seq_len(n_pcs_use), drop = FALSE]
        } else {
            idx <- sample(n_cells, min(5000, n_cells))
            dense_sub <- as.matrix(Matrix::t(centered[, idx]))
            pca_res <- stats::prcomp(dense_sub, center = FALSE,
                                     scale. = FALSE, rank. = n_pcs_use)
            pca_coords <- as.matrix(
                Matrix::t(centered) %*%
                pca_res$rotation[, seq_len(n_pcs_use), drop = FALSE]
            )
        }
    } else {
        # Standard expression-based PCA
        lib_size <- Matrix::colSums(expr_mat)
        norm_factor <- stats::median(lib_size) / lib_size

        if (nrow(expr_mat) > 2000) {
            gv <- Matrix::rowMeans(expr_mat^2) - Matrix::rowMeans(expr_mat)^2
            top_genes <- names(sort(gv, decreasing = TRUE))[seq_len(2000)]
            expr_sub <- expr_mat[top_genes, ]
        } else {
            expr_sub <- expr_mat
        }

        normed <- log1p(Matrix::t(Matrix::t(expr_sub) * norm_factor))
        gene_means_norm <- Matrix::rowMeans(normed)
        normed_centered <- normed - gene_means_norm

        n_pcs_use <- min(n_pcs,
                         min(nrow(normed_centered),
                             ncol(normed_centered)) - 1)
        if (requireNamespace("irlba", quietly = TRUE)) {
            svd_res <- irlba::irlba(Matrix::t(normed_centered),
                                    nv = n_pcs_use)
            pca_coords <- svd_res$u[, seq_len(n_pcs_use), drop = FALSE]
        } else {
            idx <- sample(n_cells, min(5000, n_cells))
            dense_sub <- as.matrix(Matrix::t(normed_centered[, idx]))
            pca_res <- stats::prcomp(dense_sub, center = FALSE,
                                     scale. = FALSE, rank. = n_pcs_use)
            pca_coords <- as.matrix(
                Matrix::t(normed_centered) %*%
                pca_res$rotation[, seq_len(n_pcs_use), drop = FALSE]
            )
        }
    }
    rownames(pca_coords) <- colnames(expr_mat)

    # Leiden clustering on SNN graph from PCA
    snn_graph <- bluster::makeSNNGraph(pca_coords, k = min(10, n_cells - 1))
    leiden_res <- igraph::cluster_leiden(
        snn_graph,
        objective_function = "modularity",
        resolution = max(1, n_clusters / 20)
    )
    assignments <- as.integer(igraph::membership(leiden_res))
    names(assignments) <- colnames(expr_mat)

    # Iteratively merge small clusters into nearest neighbor until all >= min_cells
    repeat {
        sizes <- table(assignments)
        small <- as.integer(names(sizes[sizes < min_cells]))
        if (length(small) == 0 || length(sizes) <= 1) break

        # Find the smallest cluster
        smallest_id <- as.integer(names(which.min(sizes[as.character(small)])))

        # Compute centroids for all clusters
        cluster_ids <- sort(unique(assignments))
        centers <- t(sapply(cluster_ids, function(cl) {
            colMeans(pca_coords[assignments == cl, , drop = FALSE])
        }))
        rownames(centers) <- as.character(cluster_ids)

        # Find nearest other cluster
        sc_center <- centers[as.character(smallest_id), ]
        other_ids <- setdiff(as.character(cluster_ids), as.character(smallest_id))
        other_centers <- centers[other_ids, , drop = FALSE]
        dists <- apply(other_centers, 1, function(x) sum((x - sc_center)^2))
        nearest <- as.integer(names(which.min(dists)))
        assignments[assignments == smallest_id] <- nearest
    }

    # Final check: warn if all clusters are still too small
    sizes <- table(assignments)
    if (all(sizes < min_cells)) {
        warning(
            "All clusters have fewer than min_cells (", min_cells,
            ") cells. Keeping all clusters without merging."
        )
    }

    sizes <- table(assignments)

    if (verbose) {
        message(sprintf(
            "  Final: %s microclusters, median size %s",
            length(sizes), stats::median(as.integer(sizes))
        ))
    }

    list(
        assignments = assignments,
        sizes       = sizes,
        n_clusters  = length(sizes),
        pca_coords  = pca_coords
    )
}


#' Compute cluster-inverse weights from microcluster assignments
#'
#' Each cell gets weight \code{1 / |cluster(cell)|}, so each cluster
#' contributes equally to downstream correlations regardless of size.
#' Cells in clusters smaller than \code{min_cluster_size} get weight 0.
#'
#' @param assignments Named integer vector of cluster assignments.
#' @param min_cluster_size Minimum cluster size for inclusion.
#'   Default \code{20}.
#' @return Named numeric vector of per-cell weights.
#'
#' @examples
#' sce <- sce_example()
#' mat <- SingleCellExperiment::counts(sce)
#' clusters <- create_microclusters_sparse(mat, n_clusters = 15, min_cells = 20)
#' w <- compute_simpsons_weights(clusters$assignments, min_cluster_size = 10)
#' head(w)
#'
#' @export
compute_simpsons_weights <- function(assignments, min_cluster_size = 20) {
    cluster_sizes <- table(assignments)

    w <- rep(0, length(assignments))
    names(w) <- names(assignments)

    for (cl in names(cluster_sizes)) {
        if (cluster_sizes[cl] >= min_cluster_size) {
            idx <- which(assignments == as.integer(cl))
            w[idx] <- 1.0 / cluster_sizes[cl]
        }
    }

    if (all(w == 0)) {
        warning(
            "All clusters have fewer than min_cluster_size (",
            min_cluster_size, ") cells. Using uniform weights."
        )
        w <- rep(1.0 / length(w), length(w))
        names(w) <- names(assignments)
    }

    w
}


#' Compute continuous density-based weights from PCA coordinates
#'
#' Uses kNN distance to estimate local density, then computes weights as
#' \code{median_density / local_density}. Cells in dense regions (likely
#' dominant cell types) get downweighted; cells in sparse regions get
#' upweighted. This is a continuous alternative to the discrete cluster-based
#' weights from \code{compute_simpsons_weights()}.
#'
#' @param pca_coords Matrix (cells x PCs) of PCA coordinates.
#' @param k Number of nearest neighbours for density estimation.
#'   Default \code{30}.
#' @param max_weight_quantile Quantile for capping extreme weights.
#'   Default \code{0.99}.
#' @return Named numeric vector of per-cell weights.
#'
#' @importFrom stats dist
#'
#' @keywords internal
#' @noRd
compute_density_weights <- function(pca_coords, k = 30,
                                     max_weight_quantile = 0.99) {
    n_cells <- nrow(pca_coords)
    k_use <- min(k, n_cells - 1)

    # Compute kNN distances
    if (requireNamespace("RANN", quietly = TRUE)) {
        nn <- RANN::nn2(pca_coords, k = k_use + 1)
        # Distance to k-th neighbour (skip self at index 1)
        dist_to_k <- nn$nn.dists[, k_use + 1]
    } else {
        # Fallback: brute force on subsample
        d <- as.matrix(dist(pca_coords))
        diag(d) <- Inf
        dist_to_k <- apply(d, 1, function(x) sort(x)[k_use])
    }

    # Density = 1 / distance (relative, not absolute)
    dist_to_k[dist_to_k < 1e-10] <- 1e-10
    density_est <- 1.0 / dist_to_k

    # Weights = median_density / local_density
    # Dense regions (high density) → low weight
    # Sparse regions (low density) → high weight
    med_density <- stats::median(density_est)
    w <- med_density / density_est

    # Cap extreme weights
    cap <- stats::quantile(w, max_weight_quantile)
    w[w > cap] <- cap

    # Normalize so weights sum to n_cells (preserves effective N)
    w <- w * n_cells / sum(w)
    names(w) <- rownames(pca_coords)

    w
}


#' Compute cluster-reweighted correlation matrix using FBM
#'
#' Given a normalized FBM (genes x cells) and per-cell weights,
#' computes the weighted Pearson correlation matrix as an FBM.
#' Replaces \code{extract_angles()} when Simpson's correction is active.
#'
#' @param x_mat FBM (genes x cells), already normalized.
#' @param cell_weights Numeric vector of per-cell weights
#'   (length = \code{ncol(x_mat)}).
#' @return FBM containing the weighted correlation matrix (genes x genes),
#'   with diagonal set to \code{NA}.
#'
#' @importFrom bigstatsr FBM big_apply
#'
#' @keywords internal
#' @noRd
weighted_cor_fbm <- function(x_mat, cell_weights) {
    n_genes <- nrow(x_mat)
    n_cells <- ncol(x_mat)

    # Normalize weights
    w <- cell_weights / sum(cell_weights)
    sqrt_w <- sqrt(w)

    # Step 1: Weighted means per gene
    wmeans <- bigstatsr::big_apply(
        x_mat,
        a.FUN = function(X, ind) {
            X.sub <- X[ind, , drop = FALSE]
            as.numeric(X.sub %*% w)
        },
        a.combine = "c",
        ind = bigstatsr::rows_along(x_mat),
        block.size = 500
    )

    # Step 2: Weighted SDs per gene
    wsds <- bigstatsr::big_apply(
        x_mat,
        a.FUN = function(X, ind) {
            centered <- X[ind, , drop = FALSE] - wmeans[ind]
            sqrt(as.numeric((centered^2) %*% w))
        },
        a.combine = "c",
        ind = bigstatsr::rows_along(x_mat),
        block.size = 500
    )
    wsds[wsds < .Machine$double.eps] <- NA

    # Step 3: Block-wise weighted cross-product → correlation FBM
    result <- bigstatsr::FBM(n_genes, n_genes)
    intervals <- CutBySize(n_genes, block.size = 500)

    for (j in seq_len(nrow(intervals))) {
        ind1 <- seq(intervals[j, "lower"], intervals[j, "upper"])
        block1 <- x_mat[ind1, , drop = FALSE] - wmeans[ind1]
        block1_scaled <- block1 * rep(sqrt_w, each = length(ind1))

        # Self cross-product
        cov_self <- tcrossprod(block1_scaled)
        # Convert to correlation
        sd_outer <- outer(wsds[ind1], wsds[ind1])
        cor_self <- cov_self / sd_outer
        cor_self[!is.finite(cor_self)] <- NA
        result[ind1, ind1] <- cor_self

        # Cross with remaining blocks
        if (intervals[j, "upper"] < n_genes) {
            for (k in (j + 1):nrow(intervals)) {
                ind2 <- seq(intervals[k, "lower"], intervals[k, "upper"])
                block2 <- x_mat[ind2, , drop = FALSE] - wmeans[ind2]
                block2_scaled <- block2 * rep(sqrt_w, each = length(ind2))
                cov_cross <- tcrossprod(block1_scaled, block2_scaled)
                sd_outer_cross <- outer(wsds[ind1], wsds[ind2])
                cor_cross <- cov_cross / sd_outer_cross
                cor_cross[!is.finite(cor_cross)] <- NA
                result[ind1, ind2] <- cor_cross
                result[ind2, ind1] <- t(cor_cross)
            }
        }
    }

    diag(result) <- NA
    result
}


# =============================================================================
# PART B: Local correlation preservation (cross-dataset comparison)
# =============================================================================

#' Precompute per-cluster sufficient statistics
#'
#' @param expr_mat Sparse matrix (genes x cells).
#' @param genes Character vector of gene names.
#' @param assignments Named integer vector of cluster assignments.
#' @return List with sums, sum_sq, n_cells, cluster_ids, genes.
#'
#' @importFrom Matrix rowSums
#'
#' @keywords internal
#' @noRd
precompute_simpsons_stats <- function(expr_mat, genes, assignments) {
    genes <- intersect(genes, rownames(expr_mat))
    cells <- intersect(names(assignments), colnames(expr_mat))

    expr_mat <- expr_mat[genes, cells, drop = FALSE]
    assignments <- assignments[cells]

    cluster_ids <- sort(unique(assignments))
    C <- length(cluster_ids)
    G <- length(genes)

    cluster_map <- match(assignments, cluster_ids)

    sums <- matrix(0, nrow = G, ncol = C,
                   dimnames = list(genes, cluster_ids))
    sum_sq <- matrix(0, nrow = G, ncol = C,
                     dimnames = list(genes, cluster_ids))
    n_cells <- integer(C)
    names(n_cells) <- cluster_ids

    for (ci in seq_along(cluster_ids)) {
        idx <- which(cluster_map == ci)
        n_cells[ci] <- length(idx)
        if (length(idx) == 0) next
        sub <- expr_mat[, idx, drop = FALSE]
        if (inherits(sub, "dgCMatrix") || inherits(sub, "dgTMatrix")) {
            sums[, ci] <- Matrix::rowSums(sub)
            sum_sq[, ci] <- Matrix::rowSums(sub^2)
        } else {
            sums[, ci] <- rowSums(sub)
            sum_sq[, ci] <- rowSums(sub^2)
        }
    }

    list(
        sums = sums, sum_sq = sum_sq, n_cells = n_cells,
        cluster_ids = cluster_ids, genes = genes
    )
}


#' O(1) Pearson correlation from sufficient statistics
#' @keywords internal
#' @noRd
cor_from_stats <- function(sum_x, sum_y, sum_xy, sum_x2, sum_y2, n,
                           min_n = 20) {
    if (n < min_n) return(NA_real_)
    mean_x <- sum_x / n
    mean_y <- sum_y / n
    var_x <- sum_x2 / n - mean_x^2
    var_y <- sum_y2 / n - mean_y^2
    if (var_x <= 0 || var_y <= 0) return(NA_real_)
    cov_xy <- sum_xy / n - mean_x * mean_y
    cov_xy / sqrt(var_x * var_y)
}


#' Compute per-cluster correlations for gene pairs in one batch
#'
#' @param expr_mat Sparse matrix (genes x cells).
#' @param gene_pairs Data.frame with geneA, geneB columns.
#' @param assignments Named integer vector of cluster assignments.
#' @param min_cells_cor Minimum cells per cluster for correlation.
#' @param top_quantile Quantile for peak summary.
#' @return Data.frame with geneA, geneB, peak_cor, peak_expr_level,
#'   n_valid_clusters columns.
#'
#' @importFrom stats quantile
#'
#' @keywords internal
#' @noRd
compute_local_correlations <- function(
    expr_mat,
    gene_pairs,
    assignments,
    min_cells_cor = 20,
    top_quantile = 0.95
) {
    all_genes <- unique(c(gene_pairs$geneA, gene_pairs$geneB))
    all_genes <- intersect(all_genes, rownames(expr_mat))
    cells <- intersect(names(assignments), colnames(expr_mat))
    assignments <- assignments[cells]

    cstats <- precompute_simpsons_stats(expr_mat, all_genes, assignments)
    cluster_ids <- cstats$cluster_ids

    results <- data.frame(
        geneA = gene_pairs$geneA,
        geneB = gene_pairs$geneB,
        peak_cor = NA_real_,
        peak_expr_level = NA_real_,
        n_valid_clusters = 0L,
        stringsAsFactors = FALSE
    )

    for (pi in seq_len(nrow(gene_pairs))) {
        gA <- gene_pairs$geneA[pi]
        gB <- gene_pairs$geneB[pi]

        if (!(gA %in% cstats$genes) || !(gB %in% cstats$genes)) next

        # Compute cross-products per cluster
        x <- as.numeric(expr_mat[gA, cells])
        y <- as.numeric(expr_mat[gB, cells])
        xy <- x * y
        cluster_map <- match(assignments, cluster_ids)

        cors <- numeric(0)
        expr_levels <- numeric(0)

        for (ci in seq_along(cluster_ids)) {
            n <- cstats$n_cells[ci]
            if (n < min_cells_cor) next

            idx <- which(cluster_map == ci)
            sum_xy <- sum(xy[idx])

            r <- cor_from_stats(
                sum_x = cstats$sums[gA, ci],
                sum_y = cstats$sums[gB, ci],
                sum_xy = sum_xy,
                sum_x2 = cstats$sum_sq[gA, ci],
                sum_y2 = cstats$sum_sq[gB, ci],
                n = n,
                min_n = min_cells_cor
            )

            if (!is.na(r)) {
                cors <- c(cors, abs(r))
                # Mean expression in this cluster
                mean_expr <- (cstats$sums[gA, ci] + cstats$sums[gB, ci]) /
                    (2 * n)
                expr_levels <- c(expr_levels, mean_expr)
            }
        }

        results$n_valid_clusters[pi] <- length(cors)

        if (length(cors) > 0) {
            if (length(cors) > 2) {
                peak_idx <- which(cors >= stats::quantile(cors, top_quantile))
            } else {
                peak_idx <- which.max(cors)
            }
            results$peak_cor[pi] <- mean(cors[peak_idx])
            results$peak_expr_level[pi] <- mean(expr_levels[peak_idx])
        }
    }

    results
}


#' Compare local correlation peaks across datasets
#'
#' For each gene pair, computes a Simpson's preservation score based
#' on how consistent and strong the peak local correlations are
#' across batches.
#'
#' @param local_cor_list Named list of per-batch local correlation
#'   data.frames (from \code{compute_local_correlations}).
#' @param weights Named numeric vector of batch weights.
#' @return Data.frame with geneA, geneB, simpsons_score,
#'   mean_peak_cor, sd_peak_cor columns.
#'
#' @keywords internal
#' @noRd
compare_local_correlations <- function(local_cor_list, weights) {
    # Get all gene pairs from the first batch
    template <- local_cor_list[[1]][, c("geneA", "geneB")]

    results <- data.frame(
        geneA = template$geneA,
        geneB = template$geneB,
        simpsons_score = NA_real_,
        mean_peak_cor = NA_real_,
        sd_peak_cor = NA_real_,
        stringsAsFactors = FALSE
    )

    for (pi in seq_len(nrow(template))) {
        gA <- template$geneA[pi]
        gB <- template$geneB[pi]

        peak_cors <- numeric(0)
        batch_w <- numeric(0)

        for (bn in names(local_cor_list)) {
            lcdf <- local_cor_list[[bn]]
            row_idx <- which(lcdf$geneA == gA & lcdf$geneB == gB)
            if (length(row_idx) == 1 && !is.na(lcdf$peak_cor[row_idx])) {
                peak_cors <- c(peak_cors, lcdf$peak_cor[row_idx])
                batch_w <- c(batch_w,
                             if (bn %in% names(weights)) weights[bn] else 1)
            }
        }

        if (length(peak_cors) == 0) next

        # Weighted mean and sd
        w_norm <- batch_w / sum(batch_w)
        mean_pc <- sum(w_norm * peak_cors)
        results$mean_peak_cor[pi] <- mean_pc

        if (length(peak_cors) > 1) {
            sd_pc <- sqrt(sum(w_norm * (peak_cors - mean_pc)^2))
            results$sd_peak_cor[pi] <- sd_pc
            # Score: strong AND consistent local correlations
            cv <- if (mean_pc > 0) sd_pc / mean_pc else 1
            results$simpsons_score[pi] <- mean_pc * max(0, 1 - cv)
        } else {
            results$sd_peak_cor[pi] <- 0
            results$simpsons_score[pi] <- mean_pc
        }
    }

    results
}


# =============================================================================
# PART C: Pipeline orchestration
# =============================================================================

#' Prepare Simpson's correction data for all batches
#'
#' Performs microclustering per batch and computes cluster-inverse
#' weights. Called from \code{anglemania()} before FBM conversion.
#'
#' @param matrix_list Named list of sparse matrices (genes x cells).
#' @param n_clusters Target microclusters per batch.
#' @param n_pcs PCs for clustering.
#' @param min_cells Min cells per cluster.
#' @param min_cluster_size Min cluster size for weight inclusion.
#' @param use_binary_pca Logical; cluster on the binarised detection matrix.
#' @param use_density_weights Logical; use continuous kNN density weights
#'   instead of discrete cluster-inverse weights.
#' @param density_k Number of nearest neighbours for density estimation.
#' @param use_count_split Logical; split counts into two Poisson folds so
#'   clustering and correlation use independent data.
#' @param count_split_prop Proportion of counts assigned to the clustering
#'   fold.
#' @param verbose Logical.
#' @return List with weights_per_batch and clusters_per_batch (named lists
#'   matching batch names), plus cor_matrices when count splitting is used.
#'
#' @keywords internal
#' @noRd
prepare_simpsons <- function(
    matrix_list,
    n_clusters = 15,
    n_pcs = 20,
    min_cells = 30,
    min_cluster_size = 20,
    use_binary_pca = FALSE,
    use_density_weights = FALSE,
    density_k = 30,
    use_count_split = FALSE,
    count_split_prop = 0.5,
    verbose = TRUE
) {
    weights_per_batch <- list()
    clusters_per_batch <- list()
    cor_matrices <- NULL
    if (use_count_split) {
        if (!requireNamespace("countsplit", quietly = TRUE)) {
            stop("Package 'countsplit' required for count splitting. ",
                 "Install with: install.packages('countsplit')")
        }
        cor_matrices <- list()
    }

    for (bn in names(matrix_list)) {
        if (verbose) {
            message(sprintf("  Simpson's clustering: batch '%s'...", bn))
        }

        expr_mat <- matrix_list[[bn]]

        if (use_count_split) {
            # Split counts into two independent folds
            # Fold 1 (prop) → PCA/clustering, Fold 2 (1-prop) → correlations
            split_res <- countsplit::countsplit(
                Matrix::t(expr_mat),  # countsplit expects cells x genes
                folds = 2,
                epsilon = c(count_split_prop, 1 - count_split_prop)
            )
            fold1 <- Matrix::t(split_res[[1]])  # back to genes x cells
            fold2 <- Matrix::t(split_res[[2]])
            rownames(fold1) <- rownames(fold2) <- rownames(expr_mat)
            colnames(fold1) <- colnames(fold2) <- colnames(expr_mat)

            if (verbose) {
                message(sprintf("    Count split: fold1 %d counts, fold2 %d counts",
                    sum(fold1), sum(fold2)))
            }

            # Use fold1 for clustering
            expr_mat <- fold1
            # Store fold2 for correlation computation
            cor_matrices[[bn]] <- fold2
        }

        mc <- create_microclusters_sparse(
            expr_mat = expr_mat,
            n_clusters = n_clusters,
            n_pcs = n_pcs,
            min_cells = min_cells,
            use_binary_pca = use_binary_pca,
            verbose = verbose
        )

        if (use_density_weights) {
            # Continuous density-based weights from PCA coordinates
            w <- compute_density_weights(
                pca_coords = mc$pca_coords,
                k = density_k
            )
        } else {
            # Discrete cluster-inverse weights
            w <- compute_simpsons_weights(
                assignments = mc$assignments,
                min_cluster_size = min_cluster_size
            )
        }

        weights_per_batch[[bn]] <- w
        clusters_per_batch[[bn]] <- mc$assignments
    }

    result <- list(
        weights_per_batch = weights_per_batch,
        clusters_per_batch = clusters_per_batch
    )
    if (use_count_split) {
        result$cor_matrices <- cor_matrices
    }
    result
}


#' Post-hoc Simpson's analysis on prefiltered gene pairs
#'
#' After the standard anglemania pipeline produces prefiltered gene
#' pairs, this function evaluates each pair for Simpson's paradox by
#' computing local correlations within microclusters and comparing
#' across datasets.
#'
#' @param sce SCE with anglemania results and simpsons_data in metadata.
#' @param original_matrices Named list of original sparse matrices
#'   (before FBM conversion).
#' @param min_cells_cor Minimum cells for per-cluster correlation.
#' @param top_quantile Quantile for peak summary.
#' @param verbose Logical.
#' @return Updated SCE with simpsons_score in prefiltered_df.
#'
#' @keywords internal
#' @noRd
posthoc_simpsons <- function(
    sce,
    original_matrices,
    min_cells_cor = 20,
    top_quantile = 0.95,
    verbose = TRUE
) {
    prefiltered_df <- S4Vectors::metadata(sce)$anglemania$prefiltered_df
    simpsons_data <- S4Vectors::metadata(sce)$anglemania$simpsons_data
    params <- S4Vectors::metadata(sce)$anglemania$params

    gene_pairs <- prefiltered_df[, c("geneA", "geneB")]

    if (verbose) {
        message(sprintf(
            "  Computing local correlations for %s gene pairs...",
            nrow(gene_pairs)
        ))
    }

    local_cor_list <- list()
    for (bn in names(simpsons_data$clusters_per_batch)) {
        if (verbose) {
            message(sprintf("    Batch '%s'...", bn))
        }
        local_cor_list[[bn]] <- compute_local_correlations(
            expr_mat = original_matrices[[bn]],
            gene_pairs = gene_pairs,
            assignments = simpsons_data$clusters_per_batch[[bn]],
            min_cells_cor = min_cells_cor,
            top_quantile = top_quantile
        )
    }

    # Compare across datasets
    dataset_weights <- setNames(
        params$dataset_weights$weight,
        params$dataset_weights$anglemania_batch
    )

    comparison <- compare_local_correlations(local_cor_list, dataset_weights)

    # Merge into prefiltered_df
    prefiltered_df$simpsons_score <- comparison$simpsons_score
    prefiltered_df$mean_peak_cor <- comparison$mean_peak_cor
    prefiltered_df$sd_peak_cor <- comparison$sd_peak_cor

    S4Vectors::metadata(sce)$anglemania$prefiltered_df <- prefiltered_df
    S4Vectors::metadata(sce)$anglemania$simpsons_local_cors <- local_cor_list

    sce
}
