#' Get paths to files with critical angles
#'
#' @description
#' Collects paths to the files where the critical angles were recorded to
#' during the processing of the input data by the **anglemanise** function.
#'
#' @param l_processed list. An output from the **anglemanise** function.
#' @param samp_id character vector. Samples name or names to look for.
#' @param angle_type  character vector. Specifies which type of angles to
#'   look for. Can be "sharp", "blunt", or both.
#' @return character vector. Vector of paths to the files where critical
#' angles were recorded.
#' @export get_critangs_paths
get_critangs_paths <- function(l_processed, #nolint
                               samp_id,
                               angle_type = c("sharp", "blunt")) {
  path_critangs <- paste0(
    l_processed$data_info[[samp_id]]$path_to_df_ang,
    "/",
    unlist(
      l_processed$data_info[[samp_id]][paste0("name_df_ang_", angle_type)]
    )
  )
  names(path_critangs) <- angle_type
  return(path_critangs)
}