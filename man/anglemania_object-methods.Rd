% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{show,anglemania_object-method}
\alias{show,anglemania_object-method}
\alias{matrix_list}
\alias{matrix_list,anglemania_object-method}
\alias{anglemania_object-methods,matrix_list}
\alias{matrix_list<-}
\alias{matrix_list<-,anglemania_object-method}
\alias{anglemania_object-methods,matrix_list<-}
\alias{dataset_key}
\alias{dataset_key,anglemania_object-method}
\alias{anglemania_object-methods,dataset_key}
\alias{batch_key}
\alias{batch_key,anglemania_object-method}
\alias{anglemania_object-methods,batch_key}
\alias{data_info}
\alias{data_info,anglemania_object-method}
\alias{anglemania_object-methods,data_info}
\alias{angl_weights}
\alias{angl_weights,anglemania_object-method}
\alias{anglemania_object-methods,angl_weights}
\alias{angl_weights<-}
\alias{angl_weights<-,anglemania_object-method}
\alias{anglemania_object-methods,angl_weights<-}
\alias{list_stats}
\alias{list_stats,anglemania_object-method}
\alias{anglemania_object-methods,list_stats}
\alias{list_stats<-}
\alias{list_stats<-,anglemania_object-method}
\alias{anglemania_object-methods,list_stats<-}
\alias{intersect_genes}
\alias{intersect_genes,anglemania_object-method}
\alias{anglemania_object-methods,intersect_genes}
\alias{intersect_genes<-}
\alias{intersect_genes<-,anglemania_object-method}
\alias{anglemania_object-methods,intersect_genes<-}
\alias{get_anglemania_genes}
\alias{get_anglemania_genes,anglemania_object-method}
\alias{anglemania_object-methods,get_anglemania_genes}
\alias{add_unique_batch_key}
\title{Display Summary Information for an anglemania_object}
\usage{
\S4method{show}{anglemania_object}(object)

matrix_list(object)

\S4method{matrix_list}{anglemania_object}(object)

matrix_list(object) <- value

\S4method{matrix_list}{anglemania_object}(object) <- value

dataset_key(object)

\S4method{dataset_key}{anglemania_object}(object)

batch_key(object)

\S4method{batch_key}{anglemania_object}(object)

data_info(object)

\S4method{data_info}{anglemania_object}(object)

angl_weights(object)

\S4method{angl_weights}{anglemania_object}(object)

angl_weights(object) <- value

\S4method{angl_weights}{anglemania_object}(object) <- value

list_stats(object)

\S4method{list_stats}{anglemania_object}(object)

list_stats(object) <- value

\S4method{list_stats}{anglemania_object}(object) <- value

intersect_genes(object)

\S4method{intersect_genes}{anglemania_object}(object)

intersect_genes(object) <- value

\S4method{intersect_genes}{anglemania_object}(object) <- value

get_anglemania_genes(object)

\S4method{get_anglemania_genes}{anglemania_object}(object)

add_unique_batch_key(
  object_metadata,
  dataset_key = NA_character_,
  batch_key,
  new_unique_batch_key = "batch"
)
}
\arguments{
\item{object}{An \code{anglemania_object}.}

\item{value}{A character vector of gene names.}

\item{object_metadata}{Metadata of a \code{\link[Seurat]{Seurat}} or SingleCellExperiment object.}

\item{dataset_key}{A character string specifying the column name in the
metadata that identifies the dataset. If \code{NA}, only the
\code{batch_key} is used.}

\item{batch_key}{A character string specifying the column name in the
metadata that identifies the batch.}

\item{new_unique_batch_key}{A character string for the new unique batch key
to be added to the metadata. Default is \code{"batch"}.}
}
\value{
Prints a summary to the console.

A list of \code{\link[bigstatsr]{FBM}} objects containing gene
expression matrices.

The updated \code{anglemania_object}.

A character string representing the dataset key.

A character string representing the batch key.

A data frame containing dataset and batch information.

A named numeric vector of weights.

The updated \code{anglemania_object}.

A list containing statistical matrices such as mean z-scores and SNR
z-scores

The updated \code{anglemania_object}.

A character vector of intersected gene
names from multiple Seurat or SingleCellExperiment Objects.

The updated \code{anglemania_object} object.

A character vector of integration gene names.

A \code{\link[Seurat]{Seurat}} object with an additional metadata
column containing the unique batch key.
}
\description{
This method provides a concise summary of an \code{anglemania_object}, including
dataset and batch information, the number of intersected genes, and other
relevant details.

Retrieves the list of gene expression matrices stored in the \code{anglemania_object}
object.

Assigns a new list of gene expression matrices to the \code{anglemania_object}.

Retrieves the dataset key used in the \code{anglemania_object}.

Retrieves the batch key used in the \code{anglemania_object}.

Retrieves the data frame summarizing the selected anglemania gene pairs
\code{anglemania_object}.

Retrieves the weights assigned to each dataset or batch in the \code{anglemania_object}

Assigns new weights to the datasets or batches in the \code{anglemania_object}.

Retrieves the list of statistical measures computed across datasets in the
\code{anglemania_object}.

Assigns a new list of statistical measures to the \code{anglemania_object}.

Retrieves the vector of genes that are expressed in at least the specified
number of cells across all batches.

Assigns a new vector of intersected genes to the \code{anglemania_object}.

Retrieves the list of genes selected for integration from the \code{anglemania_object}

This function adds a unique batch identifier to the metadata of a
\code{\link[Seurat]{Seurat}} object by combining specified dataset and batch
keys. This is useful for distinguishing samples during integration or
analysis.
}
\section{Functions}{
\itemize{
\item \code{show(anglemania_object)}: show anglemania_object info

\item \code{matrix_list()}: Access matrix list

\item \code{matrix_list(object) <- value}: set matrix list in anglemania_object

\item \code{dataset_key()}: Access dataset key of anglemania_object

\item \code{batch_key()}: Access batch key of anglemania_object

\item \code{data_info()}: Access info of selected gene pairs

\item \code{angl_weights()}: Access weights

\item \code{angl_weights(object) <- value}: Set weights

\item \code{list_stats()}: Access statistics of the gene-gene matrices

\item \code{list_stats(object) <- value}: Set statistics of the gene-gene matrices

\item \code{intersect_genes()}: Access the intersection of genes of all batches

\item \code{intersect_genes(object) <- value}: Set the intersection of genes of all batches

\item \code{get_anglemania_genes()}: Access the genes extracted by anglemania

\item \code{add_unique_batch_key()}: Temporarily add a unique batch key
to the dataset

}}
\examples{
sce <- sce_example()
angl <- create_anglemania_object(sce, batch_key = "batch")
angl <- anglemania(angl)
show(angl)
sce <- sce_example()
angl <- create_anglemania_object(sce, batch_key = "batch")
angl <- anglemania(angl)
str(matrix_list(angl))
sce <- sce_example()
angl <- create_anglemania_object(
  sce,
  dataset_key = "dataset",
  batch_key = "batch",
  min_cells_per_gene = 1
)
dataset_key(angl)
sce <- sce_example()
angl <- create_anglemania_object(sce, batch_key = "batch")
batch_key(angl)
sce <- sce_example()
angl <- create_anglemania_object(sce, batch_key = "batch")
batch_key(angl)
sce <- sce_example()
angl <- create_anglemania_object(sce, batch_key = "batch")
batch_key(angl)
angl_weights(angl)
sce <- sce_example()
angl <- create_anglemania_object(sce, batch_key = "batch")
angl <- anglemania(angl)
stats <- list_stats(angl)
str(stats)

sce <- sce_example()
angl <- create_anglemania_object(sce, batch_key = "batch")
intersect_genes(angl)
sce <- sce_example()
angl <- create_anglemania_object(sce, batch_key = "batch")
angl <- anglemania(angl)
# extract the genes identified by anglemania()
anglemania_genes <- get_anglemania_genes(angl)

sce <- sce_example()
angl <- create_anglemania_object(sce, batch_key = "batch")
angl <- create_anglemania_object(
  sce,
  dataset_key = "dataset",
  batch_key = "batch",
  min_cells_per_gene = 1
)
head(SingleCellExperiment::colData(sce))
}
\seealso{
\code{\link{anglemania}} \code{\link{get_list_stats}}
}
