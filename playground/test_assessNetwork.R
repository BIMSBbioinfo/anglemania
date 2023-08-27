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


build.net <- function(l_cons_nodes,
                      angle_type = "sharp") {

}