## Implementing graph search to select
## an optimum number of features for integration
#### 24.08.2023
#### ab

# Figure out the way to extract the optimum number of angles
# IE suggested to opt for a random walk;
# Perhaps a connectivity graph approach fits better to the task?
require(tidygraph)
ringgraph <- create_ring(12) #works
# Graphs are set using a tibble with 3 columns: from, to, and weights
# iterate over the recorded matrices and assemble a
# graph from conserved gene pairs. Use conservation factor as weight

##

# I must create a custom adding function that records edges
# between nodes to a file.
addmats <- function() {

}

## Approach #2 -- build a graph from recorded tables

buildnet <- function(l_processed) {

}

#### melt and filter the matrix with conserved angles
fact_sharp <- l_processed$x_sharp
fact_sharp_dt <- melt.to.df(fact_sharp)
fact_sharp_dt <- fact_sharp_dt[order(-angle)]
fringe <- 3 # select a cutting point
fact_sharp_dt <- fact_sharp_dt[angle >= fringe]
#### create a edge column for future usage
class(fact_sharp_dt) <- c("data.table", "data.frame") # dt bug here, report
fact_sharp_dt <- fact_sharp_dt[, edge := paste0(x, "_", y)]


## proceed with one sample as an example
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

samp_id <- "NB07"
paths_critangs <- get.critangs.paths(l_processed, samp_id)
## try to conjure a cmd command to filter a text file while reading
#### first, prepare a critang file in concieved format
tmp <- data.table::fread(file = paths_critangs[["sharp"]])
tmp <- tmp[, edge := paste0(x, "_", y)]
tmp <- tmp[,
           .(
             edge, x, y, angle,
             x.ref, y.ref, diff.ref, flip.ref
           )]
data.table::fwrite(tmp, file = "./critang_networktest.tsv", sep = "\t")
#### prepare a vector of conserved edges
edges_cons <- fact_sharp_dt$edge
#### filter & add an information on the node
dt_cons_samp <- data.table::fread(file = "critang_networktest.tsv", drop = 2:8)
dt_cons_samp <- dt_cons_samp[edge %chin% edges_cons]
dt_cons_samp <- dt_cons_samp[, node := samp_id]


## now time to build an algorithm;
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
          dt_cons_samp <- data.table::fread(file = critangs_paths, drop = 2:8)
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

l_cons_nodes <- assemble.cons.nodes(l_processed, c("sharp", "blunt"))



## After sleeping on the problem, can we first check relationship between
## colsum and rowsums?
names(l_cons_nodes)

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

ohm_sharp <- assemble.onehotmat(l_cons_nodes$sharp)
ohm_blunt <- assemble.onehotmat(l_cons_nodes$blunt)

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


## We can sum up both sharp and blunt angles and use the resulting
## matrix as a similarity matrix
l_conm <- list()
l_conm$sharp <- assemble.conmat(ohm_sharp)
l_conm$blunt <- assemble.conmat(ohm_blunt)
conm_summed <- (l_conm$sharp + l_conm$blunt)

# -----------------------------------------------------------------------
## what about hierarchical clustering?
hclust.conmat <- function(conm) { #nolint
  hclust(as.dist(1 - apply(conm, 2, function(x) x / diag(conm))))
}
conm_summed %>% hclust.conmat() %>% plot()

# -----------------------------------------------------------------------
## and about graph?
conm.to.tgrph <- function(conm) {
  ##
  conm <- apply(conm, 2, function(x) x / diag(conm))
  ##
  tgrph_nodes <- data.frame(name = colnames(conm))
  node_id <- seq_along(tgrph_nodes$name)
  ##
  colnames(conm) <- node_id
  rownames(conm) <- node_id
  ##
  conm[upper.tri(conm, diag = TRUE)] <- NA
  tgrph_edges <- data.table::as.data.table(conm, keep.rownames = "x")
  tgrph_edges <- data.table::melt(tgrph_edges,
                            id.vars = "x",
                            variable.name = "to",
                            value.name    = "weight")
  tgrph_edges <- tgrph_edges[!is.na(weight)]
  colnames(tgrph_edges) <- c("from", "to", "weight")
  tgrph_edges <- as.data.frame(apply(tgrph_edges, 2, as.numeric))
  ##
  tidygraph::tbl_graph(
    nodes = tgrph_nodes,
    edges = tgrph_edges,
    directed = FALSE
  )
}
## build up a graph and add some metadata to it
samps_n <- purrr::map_dbl(l_mats, ~ dim(.x)[2])
j_dmeta <- read_csv("./auxiliary/jansky_donor_meta.csv", skip = 1)
##
tgrph_summed <- conm.to.tgrph(conm_summed)
tgrph_summed <- tgrph_summed %>% 
  left_join(
    ., 
    tibble(
      name = names(samps_n),
      n = unname(samps_n)
    ),
    by = "name"
  ) %>%
  left_join(
    .,
    j_dmeta %>% dplyr::select(name = `Donor ID`, clst = `Clinical subtype`),
    by = "name"
  )

ggraph(tgrph_summed, layout = 'lgl') +
  geom_edge_link(
    aes(
      color = weight,
      width = weight
    )
   ) +
  ggraph::scale_edge_width_continuous(
    breaks = seq(0, 1, by = 0.1),
    range = c(1, 3)
  ) +
  ggraph::scale_edge_colour_gradientn(colors = viridis::inferno(n = 20)) +
  geom_node_point(aes(fill = clst), size = 8, shape = 21) +
  geom_node_label(aes(label = paste0(name, "\n", n)), repel = TRUE) +
  guides(edge_width = "none") +
  theme(legend.position = "bottom")


## Conclusion
# Graphs and hclust trees provide a great visualisation for the connectivity
# between the datasets, yet it still doesn't solve the problem with
# angle selection
save.image(file = "./snapshots/ENV_test_assessNetwork.RData")
load(file = "./snapshots/ENV_test_assessNetwork.RData")


## huinya ebanaya
### how intersection would affect gene set size?
ohm_init <- ohm_summed
tmp <- str_split_fixed(ohm_init$edge, "_", 2)
n_gns <- length(union(tmp[, 1], tmp[, 2]))
for (idx in names(diag(conm_summed))) {
    ohm_init <- ohm_init[ohm_init[[idx]] > 0, ]
    tmp <- str_split_fixed(ohm_init$edge, "_", 2)
    n_genes = length(union(tmp[, 1], tmp[, 2]))
    n_gns <- c(n_gns, n_genes)
}

### intersection cuts the data set drastically
ohm_init <- ohm_summed
##
l_ifts <- list()
tmp <- str_split_fixed(ohm_init$edge, "_", 2)
l_ifts$all <- union(tmp[, 1], tmp[, 2])
for (idx in names(diag(conm_summed))[1:4]) {
    ohm_init <- ohm_init[ohm_init[[idx]] > 0, ]
    tmp <- str_split_fixed(ohm_init$edge, "_", 2)
    l_ifts[[idx]] = union(tmp[, 1], tmp[, 2])
}

## Prepare a custom integration order
### instead of an hclust, build an integration order matrix
### from the connectivity graph
### Alternatively, build an integration order from hclust based
### on angle conservation matrix
built.int.ordr <- function(ordr_source, from = "graph") {
  if (from == "graph") {
    ##
    ordr_source %>%
    activate(edges) %>%
    data.frame() -> tmp
    purrr::map_dbl(
      1:13,
      ~ sum(tmp[tmp$from == .x | tmp$to == .x, ]$weight)
    ) -> int_order
    names(int_order) <- ordr_source %>%
      activate(nodes) %>%
      data.frame() %>%
      pull(name)
    int_order <- int_order[order(int_order, decreasing = TRUE)]
  } else if(from == "hclust") {
    int_order <- ordr_source$order
    names(int_order) <- ordr_source$labels[int_order]
  }
  ##
  for (int_step in 1:(length(int_order) - 2)) {
    ##
    if (int_step == 1) {
      binder <- rbind(
        c(-int_step, -(int_step + 1)),
        c(int_step, -(int_step + 2))
      )
    } else {
      ##
      binder <- rbind(
        binder,
        c(int_step, -(int_step + 2))
      )
    }
  }
  l_out <- list()
  l_out$int_matrix <- binder
  l_out$sample_order <- names(int_order)
  return(l_out)
}

l_int_ordr <- list()
l_int_ordr$graph <- built.int.ordr(tgrph_summed, from = "graph")
l_int_ordr$hclust <- built.int.ordr(hclust_m, from = "hclust")





## Update 26.09.2023
### Summarise the integration process and proceed to hiearachical integration

