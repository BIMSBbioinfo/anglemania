## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
####
featurise <- function(l_processed) {
  cons_nodes <- assemble.cons.nodes(
    l_processed,
    c("blunt", "sharp")
  )
  l_ufts <- purrr::map(
    cons_nodes,
    ~ count.unqiue.features(.x$conserved)
  )
  ##
  ohm <- purrr::map(
      cons_nodes,
      ~ assemble.onehotmat(.x)
    ) %>%
    purrr::reduce(., `+`)
  ##
  out <- list()
  out$unique_features <- l_ufts
  out$ohm <- ohm
  return(out)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
get.critangs.paths <- function(l_processed, #nolint
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
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
#### for every edge;node table we iterate over others and
#### append them to the graph table
# 1
## assemble a list of conserved nodes per sample
assemble.cons.nodes <- function(l_processed, angle_type, fringe = 3) { #nolint
  purrr::map(
    angle_type,
    function(angle) {
      message(paste0("Starting ", angle))
      cons_edges <- l_processed[[paste0("x_", angle)]]
      cons_edges <- melt.to.df(cons_edges)
      cons_edges <- cons_edges[order(-angle)]
      cons_edges <- cons_edges[angle >= fringe]
      cons_edges <- as.data.table(cons_edges)
      cons_edges <- cons_edges[, edge := paste0(x, "_", y)]
      cons_edges <- cons_edges[, .(edge, angle)]
      colnames(cons_edges) <- c("edge", "importance")
      ####
      nodes <- names(l_processed$data_info)
      p <- progressor(along = nodes)
      purrr::map(
        nodes,
        function(node) {
          p(message = sprintf("Processing %s", node))
          critangs_paths <- get.critangs.paths(l_processed,
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
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## Extract vectors of unique features
count.unqiue.features <- function(conserved_nodes) {
  purrr::map(
    split(conserved_nodes, conserved_nodes$importance),
    function(x) {
      str_split_fixed(x$edge, "_", 2) -> tmp
      unique(c(tmp[, 1], tmp[, 2]))
    }
  ) -> l_fts_all
  ##
  ordr <- rev(names(l_fts_all))
  res_tib <- tibble(
    importance = ordr[1],
    n_unique = length(l_fts_all[[ordr[1]]]),
    n_added = length(l_fts_all[[ordr[1]]])
  )
  for (i in 2:length(ordr)) {
     tibble(
      importance = ordr[i],
      n_unique = length(l_fts_all[[ordr[i]]]),
      n_added = length(
        intersect(
          l_fts_all[[ordr[i]]],
          l_fts_all[[ordr[i - 1]]]
        )
      )
    ) %>%
    mutate(
      n_added = n_unique - n_added
    ) -> tmp
    res_tib <- rbind(res_tib, tmp)
  }
  ##
  out <- list()
  out$feature_list <- l_fts_all
  out$features_per_ds <- res_tib
  ##
  return(out)
}

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
#### for every edge;node table we iterate over others and
#### append them to the graph table
# 1
## assemble a list of conserved nodes per sample
assemble.onehotmat <- function(cons_nodes) { #nolint
  to_keep <- str_detect(names(cons_nodes), "conserved", negate = TRUE)
  purrr::map2(
    cons_nodes[to_keep],
    names(cons_nodes[to_keep]),
    function(x, y) {
      tmp <- quote(y)
      tmp <- x[, eval(tmp) := 1]
      return(tmp)
    }
  ) %>%
    purrr::reduce(
      .,
      function(x, y) {
        merge(x, y, all = TRUE, by = "edge")
      }
    ) -> ohm
  ohm[is.na(ohm)] <- 0
  return(ohm)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## what if... bind ohms and then walk from the sample with highest number
## of critical angles down to the sample with the lowest.
# Our rationale is that the largest number of critical angles ensues
# the biggest connectivity to other samples.
assemble.conmat <- function(ohm) { #nolint
  idcs <- colnames(ohm[, 2:dim(ohm)[2]])
  ## arrange descending by the total of critical angles
  ord <- order(colSums(ohm[, ..idcs]), decreasing = TRUE)
  idcs <- idcs[ord]
  purrr::map_dfr(
    idcs,
    function(idx) {
      tmp <- ohm[, ..idcs]
      tmp <- tmp[tmp[[idx]] > 0, ]
      out <- colSums(tmp)
      return(out)
    }
  ) %>% as.matrix(.) -> m_ass
  rownames(m_ass) <- colnames(m_ass)
  return(m_ass)
}

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##
plot.unqiue.fts <- function(l_ufts, variable) {
  variable <- rlang::enquo(variable)
  purrr::map2_dfr(
    l_ufts,
    names(l_ufts),
    ~ .x$unique_features %>% mutate(angle = .y)
  ) -> t_ufts

  ggplot(
    data = t_ufts
  ) + 
  geom_line(aes(
    x = as.numeric(importance),
    y = !!variable,
    color = angle
    )
  ) + 
  scale_x_reverse()
}