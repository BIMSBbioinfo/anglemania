#'
#' @description
#' Computes the similarity of two data sets based on the conserved angles.
#' Basically, it counts the number off gene pairs present in both samples,
#' and normalizes it by the total number of the gene pairs of the samples
#' ((n_gene_pairs_samples1 + n_gene_pairs_samples2)/2).
#'
#' @param angle_list list. Output from big_anglemanise.
#'
#' @export get_integration_order
#'
get_integration_order <- function(angle_list) {
    big_dt <- rbindlist(lapply(angle_list, function(x) x$results))
    # Create a self-join to count overlaps
    overlaps <- big_dt[big_dt, on = .(gene_a, gene_b), allow.cartesian = TRUE]
    overlaps <- overlaps[dataset != i.dataset, .N, by = .(dataset, i.dataset)]
    overlap_matrix <- dcast(overlaps, dataset ~ i.dataset, value.var = "N", fill = 0)
    overlap_matrix <- overlap_matrix[, 2:ncol(overlap_matrix)]
    conm <- as.dist(1 - overlap_matrix) %>%
        hclust() %>%
        .$merge
    conm
    return(conm)
}