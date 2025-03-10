% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anglemanise_utils.R
\name{select_genes}
\alias{select_genes}
\title{Select Genes Based on Statistical Thresholds from an anglemania_object}
\usage{
select_genes(
  angl,
  zscore_mean_threshold = 2,
  zscore_sn_threshold = 2,
  max_n_genes = NULL
)
}
\arguments{
\item{angl}{An \code{anglemania_object} containing statistical
matrices such as mean z-scores and SNR z-scores.}

\item{zscore_mean_threshold}{Numeric value specifying the threshold for the
absolute mean z-score. Default is 2.}

\item{zscore_sn_threshold}{Numeric value specifying the threshold for the
SNR z-score. Default is 2.}

\item{max_n_genes}{Integer specifying the maximum number of genes to select.
If \code{NULL}, all genes that pass the thresholds are used. Default is
\code{NULL}.}
}
\value{
The input \code{anglemania_object} with the
\code{integration_genes} slot updated to include the selected genes and
their statistical information.
}
\description{
Selects genes from an \code{\link{anglemania_object-class}} based on specified
thresholds for the absolute mean z-score and signal-to-noise ratio
(SNR) z-score. It updates the \code{integration_genes} slot of the
\code{anglemania_object} with the selected genes and associated
information.
}
\details{
The function performs the following steps:
\enumerate{
\item Checks if the input object is of class \code{\link{anglemania_object-class}}.
\item If \code{max_n_genes} is not specified, it uses all genes that pass
the thresholds.
\item Identifies gene pairs where both the mean z-score and SNR z-score
exceed the specified thresholds.
\item If no gene pairs meet the criteria, it adjusts the thresholds to the
99th percentile values of the corresponding statistics and re-selects.
\item Extracts unique genes from the selected gene pairs using
\code{\link{extract_rows_for_unique_genes}}.
\item Updates the \code{integration_genes} slot of the
\code{anglemania_object}
with the selected genes and their statistics.
}
}
\examples{
sce <- sce_example()
angl <- create_anglemania_object(sce, batch_key = "batch")
angl <- anglemania(angl)
angl <- select_genes(angl,
                      zscore_mean_threshold = 2.5,
                     zscore_sn_threshold = 2.5,
                     max_n_genes = 2000)
anglemania_genes <- get_anglemania_genes(angl)
# View the selected genes and use for integration
}
\seealso{
\code{\link{extract_rows_for_unique_genes}},
\code{\link{intersect_genes}}, \code{\link{list_stats}}
}
