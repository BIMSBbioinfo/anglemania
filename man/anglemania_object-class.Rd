% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\docType{class}
\name{anglemania_object-class}
\alias{anglemania_object-class}
\alias{anglemania_object}
\title{anglemania_object - Class for Storing and Processing Gene Expression Data}
\value{
An object of class 'anglemania_object'
}
\description{
The \code{anglemania_object} class is designed to construct the correct input for the
\code{anglemania} function from a \code{\link[Seurat]{Seurat}} object and store
the results of the analysis. It encapsulates the data and metadata required
for processing gene expression data across multiple datasets and batches.
}
\section{Slots}{

\describe{
\item{\code{matrix_list}}{A list of \code{\link[bigstatsr]{FBM}} objects containing
the gene expression matrices for each batch.}

\item{\code{dataset_key}}{A character string indicating the key used to denote the
dataset in the metadata.}

\item{\code{batch_key}}{A character string indicating the key used to denote the
batch in the metadata.}

\item{\code{data_info}}{A data frame summarizing the number of samples per dataset
and their weights.}

\item{\code{weights}}{A numeric vector of weights for each dataset based on the
number of samples.}

\item{\code{list_stats}}{A list containing statistical measures computed across the
datasets.}

\item{\code{intersect_genes}}{A character vector of genes that are expressed in at
least the specified number of cells across all batches.}

\item{\code{min_cells_per_gene}}{A numeric value indicating the minimum number of
cells in which a gene must be expressed to be included in the analysis.}

\item{\code{integration_genes}}{A list containing information about integration genes
and their statistics.}

\item{\code{assay}}{A character vector specifying the Seurat assay, the counts were
extracted from.}
}}

\examples{
sce <- sce_example()
angl <- create_anglemania_object(
  sce,
  dataset_key = "dataset",
  batch_key = "batch",
  min_cells_per_gene = 1
)
angl
}
\seealso{
\code{\link{create_anglemania_object}}, \code{\link{anglemania}}
}
