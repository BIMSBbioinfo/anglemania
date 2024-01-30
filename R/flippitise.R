#' Estimate critical angles
#'
#' @description
#' A wrapper around functions that estimate_critical_angles
#' and their levels of signficance from an approximated
#' distribution of angles.
#'
#' @details
#' First, the angle distribution is approximated by quantile
#' splitting providid by the user. Then the gaussian is fit by
#' quantile transformation. Goodness of fit is tested by
#' Kolmogorov-Smirnov test. Measure of additional variance is
#' estiated fron the data by integrating the difference between
#' the fitted gaussian and theoretical distribution. Critical
#' angles are extracted from the fitted gaussian.
#'
#' @importFrom magrittr %>%
#' @param x_df_ang data.table. A long data.table with three columns:
#'   x - name of a gene in row, y - name of a gene in column,
#'   angle - an angle between x and y.
#' @param s_dims integer. Number of samples in the input gene
#'   expression matrix.
#' @param extrema double. Fraction of the angles
#'   to be cut from both sides of an approximated angle
#'   distribution.
#' @return list. Records on parameters of the gaussian,
#' goodness of fit, and critical angles.
#' @export flippitise
flippitise <- function(x_mat, #nolint
                       x_df_ang,
                       l_angles) {
  ## filter angular data.table by critical angles
  df_ang_sharp <- as.data.table(
    x_df_ang[angle <= l_angles$critical_angles[1]]
    )
  df_ang_sharp <- df_ang_sharp[, edge := paste0(x, "_", y)]
  df_ang_sharp <- df_ang_sharp[, .(edge, angle)]
  #---
  df_ang_blunt <- as.data.table(
    x_df_ang[angle >= l_angles$critical_angles[2]]
  )
  df_ang_blunt <- df_ang_blunt[, edge := paste0(x, "_", y)]
  df_ang_blunt <- df_ang_blunt[, .(edge, angle)]
  invisible(gc())
  ##
  l_flippity <- list()
  l_flippity[["sharp"]] <- df_ang_sharp
  l_flippity[["blunt"]] <- df_ang_blunt
  ##
  return(l_flippity)
}
#### Mockup --- discontinued --- requires more research
# flippitise <- function(x_mat, #nolint
#                        x_df_ang,
#                        l_angles) {
#   ## filter angular data.table by critical angles
#   df_ang_sharp <- x_df_ang[angle <= l_angles$critical_angles[1]]
#   df_ang_blunt <- x_df_ang[angle >= l_angles$critical_angles[2]]
#   invisible(gc())
#   ## prepare PCA from the input matrix
#   x_pc <- prcomp(
#     t(x_mat),
#     retx = TRUE,
#     center = FALSE,
#     scale. = FALSE,
#     rank. = 5
#   )
#   ##
#   v_featscrit <- union(
#     union(df_ang_sharp$x, df_ang_sharp$y),
#     union(df_ang_blunt$x, df_ang_blunt$y)
#   )
#   ##
#   x_mat_flippy <- rbind(
#     x_mat[v_featscrit, rownames(x_pc$x)],
#     t(x_pc$x)
#   )
#   ##
#   #tic()
#   x_mat_flippy_ang <- round(
#                             acos(
#                               1 -
#                                 parallelDist::parDist(x_mat_flippy,
#                                                       method = "cosine",
#                                                       diag = FALSE,
#                                                       upper = FALSE,
#                                                       threads = 16)
#                             ) / pi * 180,
#                             2)
#   x_mat_flippy_ang <- dist2mat(x_mat_flippy_ang, 256)
#   ##
#   r_grid <- (length(v_featscrit) + 1):(length(v_featscrit) + dim(x_pc$x)[2])
#   c_grid <- seq_along(v_featscrit)
#   x_mat_flippy_ang <- x_mat_flippy_ang[r_grid, c_grid]
#   ## doesn't matter which average to pick
#   round(colMeans(x_mat_flippy_ang), 2) -> ref
#   ## append norms to sharp and blunt angle vectors
#   ## ---------------------------------------
#   # Starting with sharp angles
#   dt_ref <- data.table::data.table(x = names(ref), x.ref = unname(ref))
#   df_ang_sharp_ref <- dt_ref[df_ang_sharp, on = "x"]
#   dt_ref <- data.table::data.table(y = names(ref), y.ref = unname(ref))
#   df_ang_sharp_ref <- dt_ref[df_ang_sharp_ref, on = "y"]
#   ## estimate parameters of the distribution of differences
#   ### For the sharp angles, we take 2 * sigma as the threshold
#   df_ang_sharp_ref <- df_ang_sharp_ref[, diff.ref := x.ref - y.ref]
#   sd_ref <- sd(df_ang_sharp_ref$diff.ref)
#   df_ang_sharp_ref <- df_ang_sharp_ref[,
#     flip.ref := ifelse(
#       abs(diff.ref) >= 2 * sd_ref,
#       sign(diff.ref),
#       0
#     )
#   ]
#   ## ---------------------------------------
#   # Proceeding with blunt angles
#   dt_ref <- data.table::data.table(x = names(ref), x.ref = unname(ref))
#   df_ang_blunt_ref <- dt_ref[df_ang_blunt, on = "x"]
#   dt_ref <- data.table::data.table(y = names(ref), y.ref = unname(ref))
#   df_ang_blunt_ref <- dt_ref[df_ang_blunt_ref, on = "y"]
#   ## Given that the distribution is centered and
#   ## bimodal, the sd would be ~= mode.
#   ## Therefore
#   df_ang_blunt_ref <- df_ang_blunt_ref[, diff.ref := x.ref - y.ref]
#   mode_ref <- sd(df_ang_blunt_ref$diff.ref)
#   df_ang_blunt_ref <- df_ang_blunt_ref[,
#     flip.ref := ifelse(
#       abs(diff.ref) >= mode_ref,
#       sign(diff.ref),
#       0
#     )
#   ]
#   ## ---------------------------------------
#   # add a node column & arrange columns in a proper way
#   df_ang_sharp_ref <- df_ang_sharp_ref[, edge := paste0(x, "_", y)]
#   df_ang_sharp_ref <- df_ang_sharp_ref[,
#                                        .(
#                                          edge, x, y, angle,
#                                          x.ref, y.ref, diff.ref, flip.ref
#                                        )]
#   df_ang_blunt_ref <- df_ang_blunt_ref[, edge := paste0(x, "_", y)]
#   df_ang_blunt_ref <- df_ang_blunt_ref[,
#                                        .(
#                                          edge, x, y, angle,
#                                          x.ref, y.ref, diff.ref, flip.ref
#                                        )]
#   l_flippity <- list()
#   l_flippity[["sharp"]] <- df_ang_sharp_ref
#   l_flippity[["blunt"]] <- df_ang_blunt_ref
#   ##
#   return(l_flippity)
# }
