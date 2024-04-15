#'
#' @description
#' Extracts the gene names from the previously generated list from big_anglemanise.
#' These genes can be used for integration later on
#'
#' @param angle_list list. Output from big_anglemanise
#' @param cutoff integer. Number of datasets where the same relationship
#'  between two genes is detected. Currently the simplest way to extract
#'  integration: setting a brute cutoff on the level of conservation
#' @return list.
#' @export big_extract_integration_features
big_extract_integration_features <- function(angle_list, cutoff) {
    union_genes <- rbindlist(lapply(angle_list, function(x) x$results)) %>%
        .[, .(signf_count = .N), by = .(gene_a, gene_b)] %>%
        .[signf_count >= cutoff] %>%
        {unique(c(.$gene_a, .$gene_b))}
}
