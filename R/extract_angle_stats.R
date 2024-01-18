#' Extract critical angles per dataset
#'
#' @description
#' Applied to the output of **anglemanise** this function will record
#' critical angles estiamted per dataset into a data.frame.
#'
#' @importFrom tidyr tibble
#' @importFrom purrr map2_dfr
#' @param l_processed list. An output from the **anglemanise** function.
#' @return data.frame. Records of critical angles per dataset.
#' @export extract_angle_stats
extract_angle_stats <- function(l_processed) {
  purrr::map2_dfr(
    l_processed[["l_angles"]],
    names(l_processed[["l_angles"]]),
    ~ cbind(
      tidyr::tibble(
        sample = .y,
        crit_sharp = .x$critical_angles[1],
        crit_blunt = .x$critical_angles[2]
      ),
      .x$statistics
    )
  )
}