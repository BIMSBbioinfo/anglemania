# ---------------------------------------------------------------------------
# Class definitions for the 'anglemania' Package
# ---------------------------------------------------------------------------
#' anglemania_object - Class for Storing and Processing Gene Expression Data
#'
#' The `anglemania_object` class is designed to construct the correct input for the
#' `anglemania` function from a \code{\link[Seurat]{Seurat}} object and store
#' the results of the analysis. It encapsulates the data and metadata required
#' for processing gene expression data across multiple datasets and batches.
#'
#' @slot matrix_list A list of \code{\link[bigstatsr]{FBM}} objects containing
#'   the gene expression matrices for each batch.
#' @slot dataset_key A character string indicating the key used to denote the
#'   dataset in the metadata.
#' @slot batch_key A character string indicating the key used to denote the
#'   batch in the metadata.
#' @slot data_info A data frame summarizing the number of samples per dataset
#'   and their weights.
#' @slot weights A numeric vector of weights for each dataset based on the
#'   number of samples.
#' @slot list_stats A list containing statistical measures computed across the
#'   datasets.
#' @slot intersect_genes A character vector of genes that are expressed in at
#'   least the specified number of cells across all batches.
#' @slot min_cells_per_gene A numeric value indicating the minimum number of
#'   cells in which a gene must be expressed to be included in the analysis.
#' @slot integration_genes A list containing information about integration genes
#'   and their statistics.
#' @slot assay A character vector specifying the Seurat assay, the counts were 
#' extracted from.
#' @return An object of class 'anglemania_object'
#' @name anglemania_object-class
#' @aliases anglemania_object
#' @docType class
#' @rdname anglemania_object-class
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(
#'   sce,
#'   dataset_key = "dataset",
#'   batch_key = "batch",
#'   min_cells_per_gene = 1
#' )
#' angl
#' @seealso \code{\link{create_anglemania_object}}, \code{\link{anglemania}}
#' @export
setClass(
  "anglemania_object",
  slots = c(
    matrix_list = "list",
    dataset_key = "character",
    batch_key = "character",
    data_info = "data.frame",
    weights = "numeric",
    list_stats = "list",
    intersect_genes = "character",
    min_cells_per_gene = "numeric",
    integration_genes = "list",
    assay = "character"
  ),
  prototype = list(
    matrix_list = list(),
    dataset_key = NA_character_,
    batch_key = NA_character_,
    data_info = data.frame(),
    weights = NA_real_,
    list_stats = list(),
    intersect_genes = NA_character_,
    min_cells_per_gene = 1,
    integration_genes = list(
      info = "data.frame",
      genes = "character"
    ),
    assay = NA_character_
  )
)

# ---------------------------------------------------------------------------
# Display Summary Information for an anglemania_object
# ---------------------------------------------------------------------------

#' Display Summary Information for an anglemania_object
#'
#' This method provides a concise summary of an \code{anglemania_object}, including
#' dataset and batch information, the number of intersected genes, and other
#' relevant details.
#'
#' @param object An \code{anglemania_object}.
#' @return Prints a summary to the console.
#' @importFrom checkmate test_string
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' angl <- anglemania(angl)
#' show(angl)
#' @describeIn anglemania_object-methods show anglemania_object info
setMethod("show", "anglemania_object", function(object) {
  cat("anglemania_object\n")
  cat("--------------\n")
  cat("Dataset key:", object@dataset_key, "\n")
  cat("Batch key:", object@batch_key, "\n")

  if (checkmate::test_string(object@dataset_key)) {
    num_datasets <- nrow(unique(object@data_info[, object@dataset_key, drop = FALSE]))
    cat("Number of datasets:", num_datasets, "\n")
    cat(
      "Datasets:",
      paste(unique(object@data_info[, object@dataset_key]), collapse = ", "),
      "\n"
    )
  } else {
    num_datasets <- 1
    cat("Number of datasets:", num_datasets, "\n")
  }

  cat("Total number of batches:", nrow(object@data_info), "\n")
  cat("Batches (showing first 5):\n")
  if (nrow(object@data_info) > 5) {
    cat(paste(object@data_info[1:5, "batch"], collapse = ", "), ", ...\n")
  } else {
    cat(paste(object@data_info[, "batch"], collapse = ", "), "\n")
  }
  cat("Number of intersected genes:", length(object@intersect_genes), "\n")
  cat("Intersected genes (showing first 10):\n")
  if (length(object@intersect_genes) > 10) {
    cat(
      paste(object@intersect_genes[1:10], collapse = ", "),
      ", ...\n"
    )
  } else {
    cat(paste(object@intersect_genes, collapse = ", "), "\n")
  }
  cat("Min cells per gene:", object@min_cells_per_gene, "\n")
})
# ---------------------------------------------------------------------------
# Accessor and Mutator Methods for the anglemania_object Class
# ---------------------------------------------------------------------------

#' Access the Matrix List from an anglemania_object
#'
#' Retrieves the list of gene expression matrices stored in the \code{anglemania_object}
#' object.
#'
#' @param object An \code{anglemania_object} object.
#' @return A list of \code{\link[bigstatsr]{FBM}} objects containing gene
#'   expression matrices.
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' angl <- anglemania(angl)
#' str(matrix_list(angl))
#' @describeIn anglemania_object-methods Access matrix list
#' @export
setGeneric(
  "matrix_list",
  function(object) standardGeneric("matrix_list")
)

#' @aliases anglemania_object-methods,matrix_list
#' @rdname anglemania_object-methods
setMethod("matrix_list", "anglemania_object", function(object) object@matrix_list)

#' Set the Matrix List in an anglemania_object
#'
#' Assigns a new list of gene expression matrices to the \code{anglemania_object}.
#'
#' @param object An \code{anglemania_object}.
#' @param value A list of \code{\link[bigstatsr]{FBM}} objects.
#' @return The updated \code{anglemania_object}.
#' @describeIn anglemania_object-methods set matrix list in anglemania_object
#' @export
setGeneric(
  "matrix_list<-",
  function(object, value) standardGeneric("matrix_list<-")
)

#' @aliases anglemania_object-methods,matrix_list<-
#' @rdname anglemania_object-methods
setReplaceMethod("matrix_list", "anglemania_object", function(object, value) {
  object@matrix_list <- value
  object
})

#' Access the Dataset Key from an anglemania_object
#'
#' Retrieves the dataset key used in the \code{anglemania_object}.
#'
#' @param object An \code{anglemania_object}.
#' @return A character string representing the dataset key.
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(
#'   sce,
#'   dataset_key = "dataset",
#'   batch_key = "batch",
#'   min_cells_per_gene = 1
#' )
#' dataset_key(angl)
#' @describeIn anglemania_object-methods Access dataset key of anglemania_object
#' @export
setGeneric(
  "dataset_key",
  function(object) standardGeneric("dataset_key")
)

#' @aliases anglemania_object-methods,dataset_key
#' @rdname anglemania_object-methods
setMethod("dataset_key", "anglemania_object", function(object) object@dataset_key)

#' Access the Batch Key from an anglemania_object
#'
#' Retrieves the batch key used in the \code{anglemania_object}.
#'
#' @param object An \code{anglemania_object}.
#' @return A character string representing the batch key.
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' batch_key(angl)
#' @describeIn anglemania_object-methods Access batch key of anglemania_object
#' @export
setGeneric(
  "batch_key",
  function(object) standardGeneric("batch_key")
)

#' @aliases anglemania_object-methods,batch_key
#' @rdname anglemania_object-methods
setMethod("batch_key", "anglemania_object", function(object) object@batch_key)

#' Access Data Information from an anglemania_object Object
#'
#' Retrieves the data frame summarizing the selected anglemania gene pairs
#' \code{anglemania_object}.
#'
#' @param object An \code{anglemania_object}.
#' @return A data frame containing dataset and batch information.
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' batch_key(angl)
#' @describeIn anglemania_object-methods Access info of selected gene pairs
#' @export
setGeneric(
  "data_info",
  function(object) standardGeneric("data_info")
)

#' @aliases anglemania_object-methods,data_info
#' @rdname anglemania_object-methods
setMethod("data_info", "anglemania_object", function(object) object@data_info)

#' Access Weights from an anglemania_object
#'
#' Retrieves the weights assigned to each dataset or batch in the \code{anglemania_object}
#'
#' @param object An \code{anglemania_object}.
#' @return A named numeric vector of weights.
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' batch_key(angl)
#' angl_weights(angl)
#' @describeIn anglemania_object-methods Access weights
#' @export
setGeneric("angl_weights", function(object) standardGeneric("angl_weights"))

#' @aliases anglemania_object-methods,angl_weights
#' @rdname anglemania_object-methods
setMethod("angl_weights", "anglemania_object", function(object) object@weights)

#' Set Weights in an anglemania_object
#'
#' Assigns new weights to the datasets or batches in the \code{anglemania_object}.
#'
#' @param object An \code{anglemania_object}.
#' @param value A named numeric vector of weights.
#' @return The updated \code{anglemania_object}.
#' @describeIn anglemania_object-methods Set weights
#' @export
setGeneric("angl_weights<-", function(object, value) standardGeneric("angl_weights<-"))

#' @aliases anglemania_object-methods,angl_weights<-
#' @rdname anglemania_object-methods
setReplaceMethod("angl_weights", "anglemania_object", function(object, value) {
  if (!is.numeric(value)) stop("weights must be numeric")
  if (is.null(names(value))) stop("weights need to be a named vector")
  # Scale the weights so that the sum of weights is 1
  message("Scaling provided weights to sum to 1...")
  value <- value / sum(value)
  object@weights <- value
  object
})

#' Access Statistical Measures from an anglemania_object
#'
#' Retrieves the list of statistical measures computed across datasets in the
#' \code{anglemania_object}.
#'
#' @param object An \code{anglemania_object}.
#' @return A list containing statistical matrices such as mean z-scores and SNR
#'   z-scores
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' angl <- anglemania(angl)
#' stats <- list_stats(angl)
#' str(stats)
#'
#' @describeIn anglemania_object-methods Access statistics of the gene-gene matrices
#' @seealso \code{\link{anglemania}} \code{\link{get_list_stats}}
#' @export
setGeneric("list_stats", function(object) standardGeneric("list_stats"))

#' @aliases anglemania_object-methods,list_stats
#' @rdname anglemania_object-methods
setMethod("list_stats", "anglemania_object", function(object) object@list_stats)

#' Set Statistical Measures in an anglemania_object
#'
#' Assigns a new list of statistical measures to the \code{anglemania_object}.
#'
#' @param object An \code{anglemania_object}.
#' @param value A list containing statistical matrices.
#' @return The updated \code{anglemania_object}.
#' @describeIn anglemania_object-methods Set statistics of the gene-gene matrices
#' @export
setGeneric("list_stats<-", function(object, value) {
  standardGeneric("list_stats<-")
})

#' @aliases anglemania_object-methods,list_stats<-
#' @rdname anglemania_object-methods
setReplaceMethod("list_stats", "anglemania_object", function(object, value) {
  if (!is.list(value)) stop("list_stats must be a list")
  object@list_stats <- value
  object
})

#' Access Intersected Genes from an anglemania_object
#'
#' Retrieves the vector of genes that are expressed in at least the specified
#' number of cells across all batches.
#'
#' @param object An \code{anglemania_object}.
#' @return A character vector of intersected gene 
#'  names from multiple Seurat or SingleCellExperiment Objects.
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' intersect_genes(angl)
#' @describeIn anglemania_object-methods 
#' Access the intersection of genes of all batches
#' @export
setGeneric(
  "intersect_genes",
  function(object) standardGeneric("intersect_genes")
)

#' @aliases anglemania_object-methods,intersect_genes
#' @rdname anglemania_object-methods
setMethod("intersect_genes", "anglemania_object", function(object) {
  object@intersect_genes
})

#' Set Intersected Genes in an anglemania_object
#'
#' Assigns a new vector of intersected genes to the \code{anglemania_object}.
#'
#' @param object An \code{anglemania_object}.
#' @param value A character vector of gene names.
#' @return The updated \code{anglemania_object} object.
#' @describeIn anglemania_object-methods 
#' Set the intersection of genes of all batches
#' @export
setGeneric("intersect_genes<-", function(object, value) {
  standardGeneric("intersect_genes<-")
})

#' @aliases anglemania_object-methods,intersect_genes<-
#' @rdname anglemania_object-methods
setReplaceMethod("intersect_genes", "anglemania_object", function(object, value) {
  object@intersect_genes <- value
  object
})

#' Extract Integration Genes from an anglemania_object
#'
#' Retrieves the list of genes selected for integration from the \code{anglemania_object}
#'
#' @param object An \code{anglemania_object}.
#' @return A character vector of integration gene names.
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' angl <- anglemania(angl)
#' # extract the genes identified by anglemania()
#' anglemania_genes <- get_anglemania_genes(angl)
#'
#' @describeIn anglemania_object-methods Access the genes extracted by anglemania
#' @export
setGeneric(
  "get_anglemania_genes",
  function(object) standardGeneric("get_anglemania_genes")
)

#' @aliases anglemania_object-methods,get_anglemania_genes
#' @rdname anglemania_object-methods
setMethod("get_anglemania_genes", "anglemania_object", function(object) {
  object@integration_genes$genes
})

#' Add a Unique Batch Key to a Seurat or SingleCellExperiment Object's Metadata
#'
#' This function adds a unique batch identifier to the metadata of a
#' \code{\link[Seurat]{Seurat}} object by combining specified dataset and batch
#' keys. This is useful for distinguishing samples during integration or
#' analysis.
#'
#' @param object_metadata Metadata of a \code{\link[Seurat]{Seurat}} or SingleCellExperiment object.
#' @param dataset_key A character string specifying the column name in the
#'   metadata that identifies the dataset. If \code{NA}, only the
#'   \code{batch_key} is used.
#' @param batch_key A character string specifying the column name in the
#'   metadata that identifies the batch.
#' @param new_unique_batch_key A character string for the new unique batch key
#'   to be added to the metadata. Default is \code{"batch"}.
#'
#' @return A \code{\link[Seurat]{Seurat}} object with an additional metadata
#'   column containing the unique batch key.
#'
#' @importFrom tidyr unite
#' @import checkmate
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' angl <- create_anglemania_object(
#'   sce,
#'   dataset_key = "dataset",
#'   batch_key = "batch",
#'   min_cells_per_gene = 1
#' )
#' head(SingleCellExperiment::colData(sce))
#' @describeIn anglemania_object-methods Temporarily add a unique batch key
#' to the dataset
#' @export
add_unique_batch_key <- function(
    object_metadata,
    dataset_key = NA_character_,
    batch_key,
    new_unique_batch_key = "batch") {

  if (checkmate::test_string(dataset_key)) {
    object_metadata <- object_metadata %>%
      tidyr::unite(
        "batch",
        dplyr::all_of(batch_key),
        sep = "_",
        remove = FALSE
      ) %>%
      tidyr::unite(
        "batch",
        dplyr::all_of(c(dataset_key, "batch")),
        sep = ":",
        remove = FALSE
      )
  } else {
    object_metadata <- object_metadata %>%
      tidyr::unite(
        "batch",
        dplyr::all_of(batch_key),
        sep = "_",
        remove = FALSE
      )
  }

  return(object_metadata)
}

# ---------------------------------------------------------------------------
# Create an anglemania_object from a Seurat or SingleCellExperiment Object
# ---------------------------------------------------------------------------

#' Create an anglemania_object from a Seurat or SingleCellExperiment Object
#'
#' Constructs an \code{\link{anglemania_object-class}} from a given
#' \code{\link[Seurat]{Seurat}} or
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' This includes extracting and processing count matrices, filtering genes
#' based on expression in a minimum number of cells, and storing results
#' along with dataset and batch information. It also calculates weights for
#' each dataset or batch based on the number of samples.
#'
#' @param object A \code{\link[Seurat]{Seurat}} or
#'    \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'    containing single-cell RNA-seq data.
#' @param dataset_key A character string indicating the column name in the
#'   object metadata that identifies the dataset to which each cell
#'   belongs. If \code{NA}, all cells are assumed to belong to the same
#'   dataset.
#' @param batch_key A character string indicating the column name(s) in the
#'   object metadata that identify the batch to which each cell belongs.
#' @param min_cells_per_gene A numeric value indicating the minimum number of
#'   cells in which a gene must be expressed to be included in the analysis.
#'   Default is \code{1}.
#' @param ... Additional arguments
#'
#' @return An \code{\link{anglemania_object-class}} containing:
#' \describe{
#'   \item{\code{matrix_list}}{A list of filtered count matrices for each unique
#'     batch.}
#'   \item{\code{dataset_key}}{The dataset key used for splitting the
#'     object.}
#'   \item{\code{batch_key}}{The batch key used for splitting the
#'     object.}
#'   \item{\code{data_info}}{A data frame summarizing the number of samples per
#'     dataset and their weights.}
#'   \item{\code{weights}}{A numeric vector of weights for each dataset or batch
#'     based on the number of samples.}
#'   \item{\code{intersect_genes}}{A character vector of genes expressed in at
#'     least the specified number of cells across all batches.}
#'   \item{\code{min_cells_per_gene}}{The minimum number of cells per gene
#'     threshold used for filtering.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Adds a unique batch key to the metadata using
#'     \code{\link{add_unique_batch_key}}.
#'   \item Extracts count matrices for each batch.
#'   \item Filters genes based on the \code{min_cells_per_gene} threshold.
#'   \item Identifies intersected genes present across all batches.
#'   \item Converts count matrices to \code{\link[bigstatsr]{FBM}} objects.
#'   \item Computes weights for each batch or dataset.
#' }
#'
#' @importFrom SeuratObject LayerData
#' @importFrom SingleCellExperiment counts
#' @importFrom tidyr unite
#' @importFrom Matrix rowSums
#' @importFrom pbapply pblapply
#' @importFrom dplyr select distinct group_by add_count mutate n_groups
#' @importFrom bigstatsr nb_cores
#' @import checkmate
#'
#' @seealso
#' \code{\link{anglemania_object-class}},
#' \code{\link{add_unique_batch_key}},
#' \code{\link{anglemania}},
#' \code{\link[bigstatsr]{FBM}}
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' angl
#' @export create_anglemania_object
#'
# Generic function
setGeneric(
  "create_anglemania_object",
  function(object, ...) standardGeneric("create_anglemania_object")
)

#' Internal helper function to get metadata from a Seurat or SCE object
#' @import checkmate
#' @importFrom SingleCellExperiment colData
#' @keywords internal
#' @noRd
.get_meta_data <- function(object) {
  if (checkmate::test_class(object, "Seurat")) {
    return(object[[]])
  } else if (checkmate::test_class(object, "SingleCellExperiment")) {
    return(SingleCellExperiment::colData(object) %>% as.data.frame())
  } else {
    stop("object must be a Seurat or SingleCellExperiment object")
  }
}

#' Internal Helper function to calculate and set weights for anglemania object
#' @importFrom dplyr select distinct group_by add_count mutate n_groups
#' @importFrom bigstatsr nb_cores
#' @keywords internal
#' @noRd
.set_weights <- function(metadata, batch_key, dataset_key = NA_character_) {
  if (checkmate::testString(dataset_key) &&
    length(unique(metadata[[dataset_key]])) > 1) {
      data_info <- metadata %>%
        dplyr::select(batch, dplyr::all_of(c(dataset_key, batch_key))) %>%
        dplyr::distinct() %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(dataset_key))) %>%
        dplyr::add_count(name = "n_samples") %>%
        dplyr::mutate(weight = 1 / n_samples / dplyr::n_groups(.)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(weight = weight / mean(weight))
    } else {
      data_info <- metadata %>%
        dplyr::select(batch, dplyr::all_of(batch_key)) %>%
        dplyr::distinct() %>%
        dplyr::mutate(weight = 1)
    }
    data_info
}


#' Internal helper function to get counts from a Seurat or SCE object
#' @import checkmate
#' @importFrom SingleCellExperiment counts
#' @importFrom SeuratObject LayerData Assays
#' @keywords internal
#' @noRd
.get_counts <- function(object, seurat_assay = "RNA", cells = NULL) {
  if (checkmate::test_class(object, "Seurat")) {
    if (seurat_assay %in% SeuratObject::Assays(object)) {
      mat <- SeuratObject::LayerData(
        object,
        assay = seurat_assay,
        layer = "counts",
        cells = cells
      )
    } else {
      stop("seurat_assay must be an assay present in the Seurat object")
    }
  } else if (checkmate::test_class(object, "SingleCellExperiment")) {
    if (is.null(cells)) {
      mat <- SingleCellExperiment::counts(object)
    } else {
      mat <- SingleCellExperiment::counts(object)[, cells, drop = FALSE]
    }
  } else {
    stop("object must be a Seurat or SingleCellExperiment object")
  }
  # check if there are any cells with 0 counts
  if (any(Matrix::colSums(mat) == 0)) {
    stop("some cells have 0 counts.\nPlease remove those cells before
     creating an anglemania object.")
  }
  return(mat)
}

#' @rdname create_anglemania_object
# Core function for creating anglemania objects
#' @param meta A data frame containing metadata for the dataset. Will be extracted from
#' the Seurat or SingleCellExperiment object or constructed from a list of objects.
#' @param matrix_list A list of count matrices for each batch.
#' @param batch_key Column name for batch information
#' @param dataset_key Column name for dataset information (optional)
#' @param min_cells_per_gene Minimum cells expressing a gene
#' @param seurat_assay Assay of the Seurat object to extract counts from
#' @param allow_missing_features Logical; allow genes missing in some batches
#' @param min_samples_per_gene Minimum samples required if
#' allow_missing_features is TRUE
#' @importFrom dplyr select distinct group_by add_count mutate n_groups
#' @importFrom Matrix rowSums
#' @importFrom pbapply pblapply
#' @importFrom bigstatsr nb_cores
#' @import checkmate
#' @return anglemania_object
#' @examples 
#' sce <- sce_example()
#' batch_key <- "batch"
#' meta <- SingleCellExperiment::colData(sce) %>% 
#'   as.data.frame() %>% 
#'   add_unique_batch_key(batch_key=batch_key)
#'
#' matrix_list <- lapply(unique(meta[[batch_key]]), function(batch) {
#'   batch_cells <- rownames(meta[meta[[batch_key]] == batch, ])
#'   batch_mat <- SingleCellExperiment::counts(sce[, batch_cells])
#'   return(batch_mat)
#' })
#'
#' angl <- .core_create_anglemania_object(
#'   meta = meta,
#'   matrix_list = matrix_list,
#'   batch_key = batch_key
#' )
#'
#' angl
#' @export
.core_create_anglemania_object <- function(
  meta,
  matrix_list,
  batch_key,
  dataset_key = NA_character_,
  min_cells_per_gene = 1,
  seurat_assay = "RNA",
  allow_missing_features = FALSE,
  min_samples_per_gene = 2
) {
  #------ CHECK INPUTS ------#
  if (checkmate::test_string(dataset_key)) {
    if (length(dataset_key) != 1) {
      stop(
        "dataset_key needs to be a character string of length 1 ",
        "corresponding to the column in the metadata of the Seurat ",
        "object that indicates which dataset the cells belong to"
      )
    } else if (!(dataset_key %in% colnames(meta))) {
      stop(
        "dataset_key needs to be a column in the metadata of the Seurat ",
        "object that indicates which dataset the cells belong to"
      )
    }
    message(
      "Using dataset_key: ",
      dataset_key
    )
  } else if (checkmate::test_scalar_na(dataset_key, null.ok = TRUE)) {
    message(
      "No dataset_key specified.\n",
      "Assuming that all samples belong to the same dataset ",
      "and are separated by batch_key: ",
      batch_key
    )
  } else {
    stop(
      "dataset_key needs to be NA/NULL or a character string of length 1 ",
      "corresponding to the column in the metadata of the Seurat ",
      "object that indicates which dataset the cells belong to"
    )
  }

  if (
    !checkmate::test_string(batch_key) ||
      length(batch_key) != 1 ||
      !(batch_key %in% colnames(meta))
  ) {
    stop(
      "batch_key needs to be a character string of length 1 ",
      "corresponding to the column in the metadata of the Seurat ",
      "object that indicates which batch the cells belong to"
    )
  }

  #------ INTERSECT GENES ------#
  # checks whether the genes should be present in all samples or just a subset of samples
  if (!allow_missing_features) {
    # Reduce to intersection of genes between batches
    message("Using the intersection of filtered genes from all batches...")
    intersect_genes <- Reduce(intersect, lapply(matrix_list, rownames))
  } else {
    message(
      "Using genes which are present in minimally ",
      min_samples_per_gene,
      " samples..."
    )
    intersect_genes <- table(unlist(lapply(matrix_list, rownames)))
    intersect_genes <- sort(rownames(intersect_genes)[
      intersect_genes >= min_samples_per_gene
    ])
  }
  message("Number of genes in intersected set: ", length(intersect_genes))
  
  #------ PADDING MATRICES & FBM CONVERSION ------#
  # Extract counts for each batch
  message("Extracting count matrices...")
  message(
    "Filtering each batch to at least ",
    min_cells_per_gene,
    " cells per gene..."
  )
  matrix_list <- pbapply::pblapply(
    matrix_list,
    function(x) {
      # Checks whether matrix contains missing features
      missing <- setdiff(intersect_genes, rownames(x))
      if (length(missing) > 0) {
        pad <- Matrix::Matrix(
          0,
          length(missing),
          ncol(x),
          dimnames = list(missing, colnames(x))
        )
        x <- rbind(x, pad)[intersect_genes, ]
        # IMPORTANT: keeps track that the missing features are ordered in the same way in all matrices
      } else {
        x <- x[match(intersect_genes, rownames(x)), ]
      }
      sparse_to_fbm(x)
    },
    cl = min(4, bigstatsr::nb_cores())
  )

  #------ DATA INFO & WEIGHTS ------#
  data_info <- .set_weights(meta, batch_key, dataset_key)
  weights <- data_info$weight
  names(weights) <- data_info$batch
  names(matrix_list) <- data_info$batch

  #------ CREATE ANGLEMANIA OBJECT ------#
  new(
    "anglemania_object",
    matrix_list = matrix_list,
    dataset_key = dataset_key,
    batch_key = batch_key,
    data_info = data_info,
    weights = weights,
    min_cells_per_gene = min_cells_per_gene,
    intersect_genes = intersect_genes,
    assay = seurat_assay
  )
}


#' @param object A Seurat object
#' @param batch_key A character string that names the column in the metadata
#' that contains the batch information
#' @param ... Additional arguments passed to .core_create_anglemania_object
#' @import checkmate
#' @importFrom Matrix rowSums
#' @export
#' @examples
#' sce <- sce_example()
#' se <- Seurat::as.Seurat(sce, data = "counts")
#' se <- SeuratObject::RenameAssays(se, "originalexp", "RNA")
#' angl <- create_anglemania_object(se, batch_key = "batch")
#' angl
#' @method create_anglemania_object Seurat

#' @rdname create_anglemania_object
#' @export
setMethod(
  "create_anglemania_object",
  "Seurat",
  function(object, batch_key, dataset_key = NA_character_, 
           min_cells_per_gene = 1, allow_missing_features = FALSE, 
           min_samples_per_gene = 2, seurat_assay = "RNA", ...) {
    
    # Metadata extraction
    meta <- .get_meta_data(object)
    meta <- add_unique_batch_key(meta, dataset_key, batch_key)
    barcodes_by_batch <- split(rownames(meta), meta$batch)

    # Matrix list creation
    matrix_list <- lapply(barcodes_by_batch, function(bc) {
      counts <- .get_counts(object, seurat_assay, cells = bc)
      counts[Matrix::rowSums(counts > 0) >= min_cells_per_gene, ]
    })

    # Call core function
    angl <- .core_create_anglemania_object(
      meta, matrix_list, batch_key, dataset_key,
      min_cells_per_gene, seurat_assay,
      allow_missing_features, min_samples_per_gene, ...
    )
    angl
  }
)

#' @rdname create_anglemania_object
#' @param object A SingleCellExperiment object
#' @param batch_key A character string that names the column in the metadata
#' that contains the batch information
#' @param ... Additional arguments passed to .core_create_anglemania_object
#' @importFrom Matrix rowSums
#' @export
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' angl
#' @method create_anglemania_object SingleCellExperiment
setMethod(
  "create_anglemania_object",
  "SingleCellExperiment",
  function(
    object,
    batch_key,
    dataset_key = NA_character_,
    min_cells_per_gene = 1,
    allow_missing_features = FALSE,
    min_samples_per_gene = 2,
    ...
  ) {
    meta <- .get_meta_data(object)
    meta <- add_unique_batch_key(meta, dataset_key, batch_key)
    barcodes_by_batch <- split(rownames(meta), meta$batch)

    matrix_list <- lapply(barcodes_by_batch, function(bc) {
      counts <- .get_counts(object, cells = bc)
      counts[Matrix::rowSums(counts > 0) >= min_cells_per_gene, ]
    })

    angl <- .core_create_anglemania_object(
      meta,
      matrix_list,
      batch_key,
      dataset_key,
      min_cells_per_gene,
      seurat_assay = "counts",
      allow_missing_features,
      min_samples_per_gene,
      ...
    )
    angl
  }
)
#' @rdname create_anglemania_object
#' @param object A list of Seurat or SingleCellExperiment objects
#' @param ... Additional arguments passed to .core_create_anglemania_object
#' @importFrom Matrix rowSums
#' @importFrom checkmate assert_list
#' @export
#' @examples
#' sce <- sce_example()
#' sce_list <- list(sce1 = sce, sce2 = sce)
#' angl <- create_anglemania_object(sce_list)
#' angl
#' se <- Seurat::as.Seurat(sce, data = "counts")
#' se <- SeuratObject::RenameAssays(se, "originalexp", "RNA")
#' se_list <- list(se1 = se, se2 = se)
#' angl <- create_anglemania_object(se_list)
#' angl
#' @method create_anglemania_object list
setMethod(
  "create_anglemania_object",
  "list",
  function(
    object,
    min_cells_per_gene = 1,
    allow_missing_features = FALSE,
    min_samples_per_gene = 2,
    seurat_assay = "RNA",
    ...
  ) {
    checkmate::assert_list(object, types = c("Seurat", "SingleCellExperiment"))

    if (is.null(names(object))) {
      names(object) <- paste0("batch", seq_along(object))
    }
    batch_key <- "batch"
    dataset_key <- NA_character_
    meta_list <- lapply(names(object), function(batch_name) {
      meta <- .get_meta_data(object[[batch_name]])
      meta$batch <- batch_name
      meta
    })
    meta <- do.call(rbind, meta_list)
    meta <- add_unique_batch_key(meta, dataset_key, batch_key)

    matrix_list <- lapply(names(object), function(batch_name) {
      obj <- object[[batch_name]]
      counts <- .get_counts(obj, seurat_assay = seurat_assay)
      counts[Matrix::rowSums(counts > 0) >= min_cells_per_gene, ]
    })

    angl <- .core_create_anglemania_object(
      meta,
      matrix_list,
      batch_key,
      dataset_key,
      min_cells_per_gene,
      seurat_assay,
      allow_missing_features,
      min_samples_per_gene,
      ...
    )
    angl
  }
)