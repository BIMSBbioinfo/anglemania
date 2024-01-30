#' Count unique features per dataset
#'
#' @description
#' Counts and records unq9iue gene features per input dataset.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @param conserved_nodes list. An output from the **assemble_cons_nodes**
#'   function.
#' @return list. First element is a list where each entry carries a
#' a vector of unique features (gene names) under sepcific conservation
#' threshold. Second element is data.frame that records the number of
#' unique features per number of datasets intersected.
#' @export count_unqiue_features
count_unqiue_features <- function(conserved_nodes) {
  purrr::map(
    split(conserved_nodes, conserved_nodes$importance),
    function(x) {
      stringr::str_split_fixed(x$edge, "_", 2) -> tmp
      unique(c(tmp[, 1], tmp[, 2]))
    }
  ) -> l_fts_all
  ##
  ##
  ordr <- rev(names(l_fts_all))
  if (length(ordr) == 1) {
    res_tib <- tidyr::tibble(
      importance = ordr,
      n_unique = length(l_fts_all[[1]]),
      n_added = length(l_fts_all[[1]]),
    )
  } else {
    res_tib <- tidyr::tibble(
      importance = ordr[1],
      n_unique = length(l_fts_all[[ordr[1]]]),
      n_added = length(l_fts_all[[ordr[1]]])
    )
    for (i in 2:length(ordr)) {
      tidyr::tibble(
        importance = ordr[i],
        n_unique = length(l_fts_all[[ordr[i]]]),
        n_added = length(
          intersect(
            l_fts_all[[ordr[i]]],
            l_fts_all[[ordr[i - 1]]]
          )
        )
      ) %>%
        dplyr::mutate(
          n_added = n_unique - n_added
        ) -> tmp
      res_tib <- rbind(res_tib, tmp)
    }
  }
  ##
  out <- list()
  out$feature_list <- l_fts_all
  out$features_per_ds <- res_tib
  ##
  return(out)
}