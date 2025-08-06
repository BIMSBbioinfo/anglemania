#' @title Adapted Reexports from bigstatsr
#' @description These are functions that were adapted from the \pkg{bigstatsr}
#' package, modified to suppress certain warnings and errors (e.g., from zero
#' variance in scaling).
#' @keywords internal
#' @name adapted_reexports
#' @rdname adapted_reexports
NULL

#' @describeIn adapted_reexports Copied version of the unexported
#' \code{bigstatsr:::CutBySize} function.
#' @importFrom bigparallelr split_len
#' @param m An integer specifying the length of the input to split
#' into intervals.
#' @param block.size An integer specifying the maximum length of each block.
#' @param nb Number of blocks. Default is \code{ceiling(m / block.size)}.
#' @return Intervals from the input length.
#' @examples
#' m = 1000
#' intervals = CutBySize(m, 100)
#' intervals
#' @export
CutBySize <- function(m, block.size, nb = ceiling(m / block.size)) {
  bigparallelr::split_len(m, nb_split = nb)
}

#' @describeIn adapted_reexports Adapted version of
#' \code{bigstatsr::big_crossprodSelf}
#' that suppresses warnings and errors related to zero scaling.
#' @importFrom bigstatsr FBM big_apply rows_along cols_along
#' big_cprodMat block_size
#' @importFrom bigparallelr seq_range
#' @inheritParams bigstatsr::big_crossprodSelf
#' @return The self crossproduct of an FBM.
#' @examples
#' mat <- matrix(rnorm(400), 20, 20)
#' X = bigstatsr::FBM(20, 20, init = mat)
#' crossp <- big_crossprodSelf_no_warning(X)
#' crossp_base <- crossprod(mat)
#' all.equal(crossp[], crossp_base)
#' @export
big_crossprodSelf_no_warning <- function(
  X,
  fun.scaling = big_scale_no_warning(center = FALSE, scale = FALSE),
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  block.size = bigstatsr::block_size(nrow(X)),
  backingfile = tempfile(tmpdir = getOption("FBM.dir"))
) {
  bigstatsr:::check_args()

  m <- length(ind.col)
  K <- bigstatsr::FBM(m, m, backingfile = backingfile)

  mu <- numeric(m)
  delta <- numeric(m)
  sums <- numeric(m)

  intervals <- CutBySize(m, block.size)

  for (j in bigstatsr::rows_along(intervals)) {
    ind1 <- bigparallelr::seq_range(intervals[j, ])
    tmp1 <- X[ind.row, ind.col[ind1]]

    ms <- fun.scaling(X, ind.row = ind.row, ind.col = ind.col[ind1])

    mu[ind1] <- ms$center
    delta[ind1] <- ms$scale
    sums[ind1] <- colSums(tmp1)

    K[ind1, ind1] <- crossprod(tmp1)

    next_lower <- intervals[j, "upper"] + 1L
    if (next_lower <= m) {
      ind2 <- next_lower:m
      K.part <- bigstatsr::big_cprodMat(
        X,
        tmp1,
        ind.row,
        ind.col[ind2],
        block.size = block.size
      )
      K[ind2, ind1] <- K.part
      K[ind1, ind2] <- t(K.part)
    }
  }

  scaleK(K, sums = sums, mu = mu, delta = delta, nrow = length(ind.row))
  structure(K, center = mu, scale = delta)
}

#' @describeIn adapted_reexports Adapted version of
#' \code{bigstatsr::big_cor} that
#' suppresses warnings and errors.
#' @inheritParams bigstatsr::big_cor
#' @return The Pearson correlation matrix from an FBM.
#' @examples
#' mat <- matrix(rnorm(400), 20, 20)
#' X = bigstatsr::FBM(20, 20, init = mat)
#' cor <- big_cor_no_warning(X)
#' cor_base <- cor(mat)
#' all.equal(cor[], cor_base)
#' @export
big_cor_no_warning <- function(
  X,
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  block.size = bigstatsr::block_size(nrow(X)),
  backingfile = tempfile(tmpdir = getOption("FBM.dir"))
) {
  cor.scaling <- function(X, ind.row, ind.col) {
    ms <- big_scale_no_warning(center = TRUE, scale = TRUE)(X, ind.row, ind.col)
    ms$scale <- ms$scale * sqrt(length(ind.row) - 1)
    ms
  }

  big_crossprodSelf_no_warning(
    X,
    fun.scaling = cor.scaling,
    ind.row = ind.row,
    ind.col = ind.col,
    block.size = block.size,
    backingfile = backingfile
  )
}

#' @describeIn adapted_reexports Adapted version of \code{bigstatsr::big_scale}
#' that suppresses warnings and errors.
#' @inheritParams bigstatsr::big_scale
#' @return A new function that returns a data.frame with two vectors,
#' \code{center} and \code{scale}, both of length \code{ind.col}.
#' @examples
#' set.seed(123)
#' mat <- matrix(rnorm(200), 20, 10)
#' X <- bigstatsr::FBM(20, 10, init = mat)
#' bs_ns <- big_scale_no_warning(center = TRUE, scale = TRUE)
#' scale_stats <- bs_ns(X)
#' scale_stats
#' bigstatsr::big_apply(X, function(X, ind) {
#'   X.sub <- X[, ind, drop = FALSE]
#'   X.sub <- t((t(X.sub) - scale_stats$center[ind]) / scale_stats$scale[ind])
#'   X[, ind] <- X.sub
#'   NULL
#' }, block.size = 3)
#' scaled_mat <- scale(mat)
#' all.equal(X[], scaled_mat, check.attributes = FALSE)
#' X[1:5, 1:5]
#' scaled_mat[1:5, 1:5]
#' @export
big_scale_no_warning <- function(center = TRUE, scale = TRUE) {
  function(
    X,
    ind.row = bigstatsr::rows_along(X),
    ind.col = bigstatsr::cols_along(X),
    ncores = 1
  ) {
    bigstatsr:::check_args()

    m <- length(ind.col)

    if (center) {
      stats <- bigstatsr::big_colstats(X, ind.row, ind.col, ncores = ncores)
      means <- stats$sum / length(ind.row)
      sds <- if (scale) sqrt(stats$var) else rep(1, m)
    } else {
      means <- rep(0, m)
      sds <- rep(1, m)
    }

    data.frame(center = means, scale = sds)
  }
}
