#' Extract integration features (gene names)
#'
#' @description
#' Number of gene pairs with conserved angles diminish along the number
#' of intersected dataset, seemingly resembling a Poisson distribution.
#' Here, we use the Poisson distribution to provide some sort of background
#' to the selection of the cutoff for the number of the intersected samples.
#'
#' @import data.table
#' @import tidyr
#' @import dplyr
#' @importFrom magrittr %>%
#' @param l_cumfact list. First element of the output from the
#'   **anglemanise** function.
#' @return list. Two elements contain character vector with gene names
#' corresponding to sharp and blutn angles.
#' @export
extract.integration.features <- function(l_cumfact, #nolint
                                         cutoff = NULL) {
  purrr::map(
    l_cumfact[1:2],
    function(cumfact_mat) {
      cumfact_dt <- data.table::as.data.table(cumfact_mat, keep.rownames = "x")
      cumfact_dt_melt <- data.table::melt(
        cumfact_dt,
        id.vars = "x",
        variable.name = "y",
        value.name    = "angle.fact"
      )
      ## Clean the table
      cumfact_dt_melt <- cumfact_dt_melt[!is.na(angle.fact) & (x != y)]
      rm(cumfact_dt); gc()
      ## this distribution likely follows poisson ??? doubts
      if (is.null(cutoff)) {
        cumfact_dt_melt %>%
          janitor::tabyl(angle.fact) %>%
          tidyr::as_tibble() -> fact_stats
        glm(
          n ~ angle.fact,
          data = fact_stats,
          family = quasi(variance = "mu", link = "log")
        ) -> glmodel.poisson
        fact_stats %>%
          dplyr::mutate(
            n_pred_poisson = predict(
              glmodel.poisson,
              newdata = tidyr::tibble(angle.fact = fact_stats$angle.fact)
            ) %>%
              exp() %>%
              unname() %>%
              as.integer()
          ) -> fact_stats_pred
        fact_stats_pred %>%
          dplyr::transmute(
            diff = log(n/n_pred_poisson),
            diff = ifelse(diff == Inf, 0, diff)
          ) %>%
          dplyr::pull(diff) -> v_diff
        which(v_diff == max(v_diff)) - 1 -> cutoff
      } # first -1 to get to 0-based counting;
      ##
      cumfact_dt_melt_filt <- cumfact_dt_melt[angle.fact >= cutoff]
      rm(cumfact_dt_melt); gc()
      features <- intersect(cumfact_dt_melt_filt$x, cumfact_dt_melt_filt$y)
      ##
      return(features)
    }
  ) -> l_features
  ##
  names(l_features) <- names(l_cumfact[1:2])
  return(l_features)
}