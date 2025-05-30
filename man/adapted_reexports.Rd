% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexports_bigstatsr.R
\name{adapted_reexports}
\alias{adapted_reexports}
\alias{CutBySize}
\alias{big_crossprodSelf_no_warning}
\alias{big_cor_no_warning}
\alias{big_scale_no_warning}
\title{Adapted Reexports from bigstatsr}
\usage{
CutBySize(m, block.size, nb = ceiling(m/block.size))

big_crossprodSelf_no_warning(
  X,
  fun.scaling = big_scale_no_warning(center = FALSE, scale = FALSE),
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  block.size = bigstatsr::block_size(nrow(X)),
  backingfile = tempfile(tmpdir = getOption("FBM.dir"))
)

big_cor_no_warning(
  X,
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  block.size = bigstatsr::block_size(nrow(X)),
  backingfile = tempfile(tmpdir = getOption("FBM.dir"))
)

big_scale_no_warning(center = TRUE, scale = TRUE)
}
\arguments{
\item{m}{An integer specifying the length of the input to split into intervals.}

\item{block.size}{An integer specifying the maximum length of each block.}

\item{nb}{Number of blocks. Default is \code{ceiling(m / block.size)}.}

\item{X}{An object of class \link[bigstatsr:FBM-class]{FBM}.}

\item{fun.scaling}{A function with parameters \code{X}, \code{ind.row} and \code{ind.col},
and that returns a data.frame with \verb{$center} and \verb{$scale} for the columns
corresponding to \code{ind.col}, to scale each of their elements such as followed:
\deqn{\frac{X_{i,j} - center_j}{scale_j}.} Default doesn't use any scaling.
You can also provide your own \code{center} and \code{scale} by using \code{\link[bigstatsr:as_scaling_fun]{as_scaling_fun()}}.}

\item{ind.row}{An optional vector of the row indices that are used.
If not specified, all rows are used. \strong{Don't use negative indices.}}

\item{ind.col}{An optional vector of the column indices that are used.
If not specified, all columns are used. \strong{Don't use negative indices.}}

\item{backingfile}{Path to the file storing the FBM data on disk.
\strong{An extension ".bk" will be automatically added.} Default stores in the
temporary directory, which you can change using global option "FBM.dir".}

\item{center}{A logical value: whether to return means or 0s.}

\item{scale}{A logical value: whether to return standard deviations or 1s.
\strong{You can't use scale without using center.}}
}
\value{
Intervals from the input length.

The self crossproduct of an FBM.

The Pearson correlation matrix from an FBM.

A new function that returns a data.frame with two vectors, \code{center}
and \code{scale}, both of length \code{ind.col}.
}
\description{
These are functions that were adapted from the \pkg{bigstatsr}
package, modified to suppress certain warnings and errors (e.g., from zero
variance in scaling).
}
\section{Functions}{
\itemize{
\item \code{CutBySize()}: Copied version of the unexported
\code{bigstatsr:::CutBySize} function.

\item \code{big_crossprodSelf_no_warning()}: Adapted version of \code{bigstatsr::big_crossprodSelf}
that suppresses warnings and errors related to zero scaling.

\item \code{big_cor_no_warning()}: Adapted version of \code{bigstatsr::big_cor} that
suppresses warnings and errors.

\item \code{big_scale_no_warning()}: Adapted version of \code{bigstatsr::big_scale}
that suppresses warnings and errors.

}}
\examples{
m = 1000
intervals = CutBySize(m, 100)
intervals
mat <- matrix(rnorm(400), 20, 20)
X = bigstatsr::FBM(20, 20, init = mat)
crossp <- big_crossprodSelf_no_warning(X)
crossp_base <- crossprod(mat)
all.equal(crossp[], crossp_base)
mat <- matrix(rnorm(400), 20, 20)
X = bigstatsr::FBM(20, 20, init = mat)
cor <- big_cor_no_warning(X)
cor_base <- cor(mat)
all.equal(cor[], cor_base)
set.seed(123)
mat <- matrix(rnorm(200), 20, 10)
X <- bigstatsr::FBM(20, 10, init = mat)
bs_ns <- big_scale_no_warning(center = TRUE, scale = TRUE)
scale_stats <- bs_ns(X)
scale_stats
bigstatsr::big_apply(X, function(X, ind) {
  X.sub <- X[, ind, drop = FALSE]
  X.sub <- t((t(X.sub) - scale_stats$center[ind]) / scale_stats$scale[ind])
  X[, ind] <- X.sub
  NULL
}, block.size = 3)
scaled_mat <- scale(mat)
all.equal(X[], scaled_mat, check.attributes = FALSE)
X[1:5, 1:5]
scaled_mat[1:5, 1:5]
}
\keyword{internal}
