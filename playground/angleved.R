suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
  library(ggforce)
  library(matrixStats)
})

plot_vectors = function(
    dv = data.frame(X = c(2,1), Y = c(1,2)),
    radius    = 2,
    linewidth = 1.5,
    slope = FALSE,
    plot = TRUE,
    title = ""
){
  if(class(dv) != "data.frame")
    dv = data.frame(dv) 
  
  dv = dv %>%
    magrittr::set_colnames(c("X", "Y")) %>%
    mutate(X0 = 0) %>%
    mutate(Y0 = 0) %>%
    mutate(col = factor(1:n()))
  
  gg = dv %>%
    ggplot(aes(x = X0, y = Y0, xend = X, yend = Y, color=col)) +
    geom_segment(arrow = arrow(length = unit(0.03, "npc")), linewidth=linewidth) + 
    theme_bw() + 
    xlim(c(-radius,radius)) +
    ylim(c(-radius,radius)) +
    geom_hline(yintercept = 0, linewidth=linewidth-0.75, color="gray") +
    geom_vline(xintercept = 0, linewidth=linewidth-0.75, color="gray") + 
    theme(legend.position="none") +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), color="black")

  if(slope)
    gg = gg + geom_abline(slope = 1)
  
  if(nchar(title) > 0)
    gg = gg + ggtitle(title)
  
  if(plot){
    print(gg)
  }else{
   return(gg) 
  }
  
}

rotate_vector = function(
    v = c(1,2),
    theta = 45
){
  rotation_matrix = matrix(c(cos(pi*theta/180), sin(pi*theta/180), -sin(pi*theta/180), cos(pi*theta/360)), nrow=2)
  rotation_matrix %*% v
}

rotate_vectors = function(
    v, 
    thetas = 45
){
  lv = lapply(thetas, function(x) rotate_vector(v = v, theta=x))
  lv = lapply(lv, t)
  dv = do.call(rbind, lv) %>%
      as.data.frame() %>%
      magrittr::set_colnames(c("X","Y"))
  dv
}

#
## prepare graph solutions
## make it into a package
## play with norm vs scaled vectors
#


## simulate a random matrix (negative correlations present)
library(tidyverse)
melt.to.df <- function(x_mat_ang) { #nolint
  ##
  x_mat_ang[upper.tri(x_mat_ang, diag = TRUE)] <- NA
  x_df_ang <- data.table::as.data.table(x_mat_ang, keep.rownames = "x")
  x_df_ang <- data.table::melt(x_df_ang,
                               id.vars = "x",
                               variable.name = "y",
                               value.name    = "angle")
  x_df_ang <- na.omit(x_df_ang)
  x_df_ang <- x_df_ang[, y := as.character(y)]
  class(x_df_ang) <- "data.frame"
  return(x_df_ang)
}
# Number of features
n <- 10
m_simexp <- purrr::map_dfc(
  1:n,
  ~ tibble(!!sym(paste0("gene_", .x)) := rnegbin(15, 12, 0.4))
) %>% 
as.matrix() %>% t() %>% `colnames<-`(paste0("cell_", 1:15))

get.ang.dt <- function(m_exp) {
    m_ang <- acos(
    1 -
        parallelDist::parDist(m_exp,
                            method = "cosine",
                            diag = FALSE,
                            upper = FALSE,
                            threads = 16)
    ) / pi * 180
    melt.to.df(x_mat_ang = as.matrix(m_ang))
}

l_ass <- list()
l_ass$raw  <- get.ang.dt(m_simexp)
l_ass$norm <- get.ang.dt(t(apply(m_simexp, 1, function(x) x / sum(abs(x)))))
l_ass$scaled <- get.ang.dt(t(apply(m_simexp, 1, function(x) (x - mean(x)) / sd(x))))
purrr::map(l_ass,
           ~ ggplot(data = .x, aes(x = angle)) + 
           geom_histogram() +
           scale_x_continuous(breaks = seq(0, 180, by = 15), 
                              limits = c(0, 180)
                              )
            )

## makes sense
## how does the angles change then?
inner_join(
    l_ass$raw,
    l_ass$scaled,
    by = c("x", "y")
) %>% `colnames<-`(c("x", "y", "raw", "scaled")) %>%
ggplot(data = ., aes(x = raw, y = scaled)) + 
geom_point() + 
xlim(c(0, 90)) + ylim(c(0, 180)) +
geom_abline(intercept = 0, slope = 1) +
theme(aspect.ratio = 1)
save.image(file = "./snapshots/ENV_angleved.RData")


