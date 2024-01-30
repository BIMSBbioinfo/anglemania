#' Assebmle gene pairs with conserved angles
#'
#' @description
#' Extract unique gene features per dataset with significant
#' sharp/blunt angles.
#'
#' @importFrom data.table as.data.table fread
#' @import progressr
#' @importFrom purrr map
#' @param l_processed list. An output from the **anglemanise** function.
#' @param angle_type character vector. Specifies which type of angles to
#'   assemble. Can be "sharp", "blunt", or both.
#' @param fringe integer. Minimal number of samples where the relationship
#'   between genes should be preserved, otherwise the gene pair is removed.
#'   Default is 3.
#' @return list. Each element contains a list with entries recording
#' significant angles per sample and an entry recording conserved angles
#' and their conservation level.
#' @export assemble_cons_nodes
assemble_cons_nodes <- function(l_processed, angle_type, fringe = 3) { #nolint
  purrr::map(
    angle_type,
    function(angle) {
      message(paste0("Starting ", angle))
      cons_edges <- l_processed[[paste0("x_", angle)]]
      cons_edges <- melt_to_df(cons_edges)
      cons_edges <- cons_edges[order(-angle)]
      cons_edges <- cons_edges[angle >= fringe]
      cons_edges <- data.table::as.data.table(cons_edges)
      cons_edges <- cons_edges[, edge := paste0(x, "_", y)]
      cons_edges <- cons_edges[, .(edge, angle)]
      colnames(cons_edges) <- c("edge", "importance")
      ####
      nodes <- names(l_processed$data_info)
      p <- progressr::progressor(along = nodes)
      purrr::map(
        nodes,
        function(node) {
          p(message = sprintf("Processing %s", node))
          critangs_paths <- get_critangs_paths(l_processed,
                                                          node,
                                                          angle)
          dt_cons_samp <- data.table::fread(file = critangs_paths, drop = 2)
          dt_cons_samp <- dt_cons_samp[edge %chin% cons_edges$edge]
        }
      ) -> l_nodes_samp
      names(l_nodes_samp) <- nodes
      l_nodes_samp <- c(l_nodes_samp, conserved = list(cons_edges))
    }
  ) -> l_nodes
  names(l_nodes) <- angle_type
  return(l_nodes)
}