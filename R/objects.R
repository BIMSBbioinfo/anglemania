# ---------------------------------------------------------------------------
# Class definitions for the 'anglemania' Package
# ---------------------------------------------------------------------------
#' anglemaniaObject - Class for Storing and Processing Gene Expression Data
#'
#' The `anglemaniaObject` class is designed to construct the correct input for the
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
#'
#' @name anglemaniaObject-class
#' @rdname anglemaniaObject-class
#' @examples
#' load(system.file(
#'   "extdata",
#'   "seurat_splatter_sim.RData",
#'   package = "anglemania"
#' ))
#'
#' se[[]]$Dataset <- rep(c("A", "B"), each = ncol(se) / 2)
#' anglemania_object <- create_anglemaniaObject(
#'   se,
#'   dataset_key = "Dataset",
#'   batch_key = "Batch",
#'   min_cells_per_gene = 1
#' )
#' anglemania_object
#' @seealso \code{\link{create_anglemaniaObject}}, \code{\link{anglemania}}
#' @exportClass anglemaniaObject
setClass(
  "anglemaniaObject",
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
# Display Summary Information for an anglemaniaObject
# ---------------------------------------------------------------------------

#' Display Summary Information for an anglemaniaObject
#'
#' This method provides a concise summary of an \code{anglemaniaObject}, including
#' dataset and batch information, the number of intersected genes, and other
#' relevant details.
#'
#' @param object An \code{anglemaniaObject}.
#' @return Prints a summary to the console.
#' @importFrom checkmate testString
#' @examples
#' load(system.file(
#'   "extdata",
#'   "seurat_splatter_sim.RData",
#'   package = "anglemania"
#' ))
#'
#' anglemania_object <- create_anglemaniaObject(
#'   se,
#'   batch_key = "Batch",
#'   min_cells_per_gene = 1
#' )
#' show(anglemania_object)
#' @describeIn anglemaniaObject-methods show anglemaniaObject info
setMethod("show", "anglemaniaObject", function(object) {
  cat("anglemaniaObject\n")
  cat("--------------\n")
  cat("Dataset key:", object@dataset_key, "\n")
  cat("Batch key:", object@batch_key, "\n")

  if (checkmate::testString(object@dataset_key)) {
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
# Accessor and Mutator Methods for the anglemaniaObject Class
# ---------------------------------------------------------------------------

#' Access the Matrix List from an anglemaniaObject
#'
#' Retrieves the list of gene expression matrices stored in the \code{anglemaniaObject}
#' object.
#'
#' @param object An \code{anglemaniaObject} object.
#' @return A list of \code{\link[bigstatsr]{FBM}} objects containing gene
#'   expression matrices.
#' @examples
#' load(system.file(
#'   "extdata",
#'   "seurat_splatter_sim.RData",
#'   package = "anglemania"
#' ))
#'
#' anglemania_object <- create_anglemaniaObject(
#'   se,
#'   batch_key = "Batch",
#'   min_cells_per_gene = 1
#' )
#' str(matrix_list(anglemania_object))
#' @describeIn anglemaniaObject-methods Access matrix list
#' @export
setGeneric(
  "matrix_list",
  function(object) standardGeneric("matrix_list")
)
setMethod("matrix_list", "anglemaniaObject", function(object) object@matrix_list)

#' Set the Matrix List in an anglemaniaObject
#'
#' Assigns a new list of gene expression matrices to the \code{anglemaniaObject}.
#'
#' @param object An \code{anglemaniaObject}.
#' @param value A list of \code{\link[bigstatsr]{FBM}} objects.
#' @return The updated \code{anglemaniaObject}.
#' @describeIn anglemaniaObject-methods set matrix list in anglemaniaObject
#' @keywords internal
setGeneric(
  "matrix_list<-",
  function(object, value) standardGeneric("matrix_list<-")
)
setReplaceMethod("matrix_list", "anglemaniaObject", function(object, value) {
  object@matrix_list <- value
  object
})

#' Access the Dataset Key from an anglemaniaObject
#'
#' Retrieves the dataset key used in the \code{anglemaniaObject}.
#'
#' @param object An \code{anglemaniaObject}.
#' @return A character string representing the dataset key.
#' @examples
#' load(system.file(
#'   "extdata",
#'   "seurat_splatter_sim.RData",
#'   package = "anglemania"
#' ))
#'
#' anglemania_object <- create_anglemaniaObject(
#'   se,
#'   batch_key = "Batch",
#'   min_cells_per_gene = 1
#' )
#' dataset_key(anglemania_object)
#' @describeIn anglemaniaObject-methods Access dataset key of anglemaniaObject
#' @export
setGeneric(
  "dataset_key",
  function(object) standardGeneric("dataset_key")
)
setMethod("dataset_key", "anglemaniaObject", function(object) object@dataset_key)

#' Access the Batch Key from an anglemaniaObject
#'
#' Retrieves the batch key used in the \code{anglemaniaObject}.
#'
#' @param object An \code{anglemaniaObject}.
#' @return A character string representing the batch key.
#' @examples
#' load(system.file(
#'   "extdata",
#'   "seurat_splatter_sim.RData",
#'   package = "anglemania"
#' ))
#'
#' anglemania_object <- create_anglemaniaObject(
#'   se,
#'   batch_key = "Batch",
#'   min_cells_per_gene = 1
#' )
#' batch_key(anglemania_object)
#' @describeIn anglemaniaObject-methods Access batch key of anglemaniaObject
#' @export
setGeneric(
  "batch_key",
  function(object) standardGeneric("batch_key")
)
setMethod("batch_key", "anglemaniaObject", function(object) object@batch_key)

#' Access Data Information from an anglemaniaObject Object
#'
#' Retrieves the data frame summarizing the selected anglemania gene pairs
#' \code{anglemaniaObject}.
#'
#' @param object An \code{anglemaniaObject}.
#' @return A data frame containing dataset and batch information.
#' @examples
#' load(system.file(
#'   "extdata",
#'   "seurat_splatter_sim.RData",
#'   package = "anglemania"
#' ))
#'
#' anglemania_object <- create_anglemaniaObject(
#'   se,
#'   batch_key = "Batch",
#'   min_cells_per_gene = 1
#' )
#' data_info(anglemania_object)
#' @describeIn anglemaniaObject-methods Access info of selected gene pairs
#' @export
setGeneric(
  "data_info",
  function(object) standardGeneric("data_info")
)
setMethod("data_info", "anglemaniaObject", function(object) object@data_info)

#' Access Weights from an anglemaniaObject
#'
#' Retrieves the weights assigned to each dataset or batch in the \code{anglemaniaObject}
#'
#' @param object An \code{anglemaniaObject}.
#' @return A named numeric vector of weights.
#' @examples
#' load(system.file(
#'   "extdata",
#'   "seurat_splatter_sim.RData",
#'   package = "anglemania"
#' ))
#'
#' anglemania_object <- create_anglemaniaObject(
#'   se,
#'   batch_key = "Batch",
#'   min_cells_per_gene = 1
#' )
#' angl_weights(anglemania_object)
#' @describeIn anglemaniaObject-methods Access weights
#' @export
setGeneric("angl_weights", function(object) standardGeneric("angl_weights"))
setMethod("angl_weights", "anglemaniaObject", function(object) object@weights)

#' Set Weights in an anglemaniaObject
#'
#' Assigns new weights to the datasets or batches in the \code{anglemaniaObject}.
#'
#' @param object An \code{anglemaniaObject}.
#' @param value A named numeric vector of weights.
#' @return The updated \code{anglemaniaObject}.
#' @describeIn anglemaniaObject-methods Set weights
#' @keywords internal
setGeneric("angl_weights<-", function(object, value) standardGeneric("angl_weights<-"))
setReplaceMethod("angl_weights", "anglemaniaObject", function(object, value) {
  if (!is.numeric(value)) stop("weights must be numeric")
  if (is.null(names(value))) stop("weights need to be a named vector")
  # Scale the weights so that the sum of weights is 1
  message("Scaling provided weights to sum to 1...")
  value <- value / sum(value)
  object@weights <- value
  object
})

#' Access Statistical Measures from an anglemaniaObject
#'
#' Retrieves the list of statistical measures computed across datasets in the
#' \code{anglemaniaObject}.
#'
#' @param object An \code{anglemaniaObject}.
#' @return A list containing statistical matrices such as mean z-scores and SNR
#'   z-scores
#' @examples
#'
#' # list_stats extracts the statistical measures from the anglemaniaObject
#' # after running anglemania()
#' stats <- list_stats(anglemania_object)
#'
#' @describeIn anglemaniaObject-methods Access statistics of the gene-gene matrices
#' @seealso \code{\link{anglemania}} \code{\link{get_list_stats}}
#' @export
setGeneric("list_stats", function(object) standardGeneric("list_stats"))
setMethod("list_stats", "anglemaniaObject", function(object) object@list_stats)

#' Set Statistical Measures in an anglemaniaObject
#'
#' Assigns a new list of statistical measures to the \code{anglemaniaObject}.
#'
#' @param object An \code{anglemaniaObject}.
#' @param value A list containing statistical matrices.
#' @return The updated \code{anglemaniaObject}.
#' @describeIn anglemaniaObject-methods Set statistics of the gene-gene matrices
#' @keywords internal
setGeneric("list_stats<-", function(object, value) {
  standardGeneric("list_stats<-")
})
setReplaceMethod("list_stats", "anglemaniaObject", function(object, value) {
  if (!is.list(value)) stop("list_stats must be a list")
  object@list_stats <- value
  object
})

#' Access Intersected Genes from an anglemaniaObject
#'
#' Retrieves the vector of genes that are expressed in at least the specified
#' number of cells across all batches.
#'
#' @param object An \code{anglemaniaObject}.
#' @return A character vector of intersected gene 
#'  names from multiple Seurat objects.
#' @examples
#' load(system.file(
#'   "extdata",
#'   "seurat_splatter_sim.RData",
#'   package = "anglemania"
#' ))
#'
#' anglemania_object <- create_anglemaniaObject(
#'   se,
#'   batch_key = "Batch",
#'   min_cells_per_gene = 1
#' )
#' intersect_genes(anglemania_object)
#' @describeIn anglemaniaObject-methods 
#' Access the intersection of genes of all batches
#' @export
setGeneric(
  "intersect_genes",
  function(object) standardGeneric("intersect_genes")
)
setMethod("intersect_genes", "anglemaniaObject", function(object) {
  object@intersect_genes
})

#' Set Intersected Genes in an anglemaniaObject
#'
#' Assigns a new vector of intersected genes to the \code{anglemaniaObject}.
#'
#' @param object An \code{anglemaniaObject}.
#' @param value A character vector of gene names.
#' @return The updated \code{anglemaniaObject} object.
#' @describeIn anglemaniaObject-methods 
#' Set the intersection of genes of all batches
#' @keywords internal
setGeneric("intersect_genes<-", function(object, value) {
  standardGeneric("intersect_genes<-")
})
setReplaceMethod("intersect_genes", "anglemaniaObject", function(object, value) {
  object@intersect_genes <- value
  object
})

#' Extract Integration Genes from an anglemaniaObject
#'
#' Retrieves the list of genes selected for integration from the \code{anglemaniaObject}
#'
#' @param object An \code{anglemaniaObject}.
#' @return A character vector of integration gene names.
#' @examples
#'
#' # extract the genes identified by anglemania()
#' anglemania_genes <- get_anglemania_genes(anglemania_object)
#'
#' @describeIn anglemaniaObject-methods Access the genes extracted by anglemania
#' @export
setGeneric(
  "get_anglemania_genes",
  function(object) standardGeneric("get_anglemania_genes")
)
setMethod("get_anglemania_genes", "anglemaniaObject", function(object) {
  object@integration_genes$genes
})

#' Add a Unique Batch Key to a Seurat Object's Metadata
#'
#' This function adds a unique batch identifier to the metadata of a
#' \code{\link[Seurat]{Seurat}} object by combining specified dataset and batch
#' keys. This is useful for distinguishing samples during integration or
#' analysis.
#'
#' @param seurat_object A \code{\link[Seurat]{Seurat}} object.
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
#' @examples 
#' load(system.file(
#'  "extdata",
#'  "seurat_splatter_sim.RData",
#'  package = "anglemania"))
#' 
#' se[[]]$Dataset <- rep(c("A", "B"), each = ncol(se)/2)
#' seurat_object <- add_unique_batch_key(
#'   seurat_object = se,
#'   dataset_key = "Dataset",
#'   batch_key = "Batch",
#'   new_unique_batch_key = "batch" 
#'   )
#' head(seurat_object[[]])
#' @describeIn anglemaniaObject-methods Temporarily add a unique batch key to the dataset
#' @export
add_unique_batch_key <- function(
    seurat_object,
    dataset_key = NA_character_,
    batch_key,
    new_unique_batch_key = "batch") {
  if (checkmate::testString(dataset_key)) {
    meta <- seurat_object[[]] %>%
      tidyr::unite(
        "batch",
        all_of(batch_key),
        sep = "_",
        remove = FALSE
      ) %>%
      tidyr::unite(
        "batch",
        all_of(c(dataset_key, "batch")),
        sep = ":",
        remove = FALSE
      )
  } else {
    meta <- seurat_object[[]] %>%
      tidyr::unite(
        "batch",
        all_of(batch_key),
        sep = "_",
        remove = FALSE
      )
  }
  seurat_object[[]] <- meta
  return(seurat_object)
}

# ---------------------------------------------------------------------------
# Create an anglemaniaObject from a Seurat Object
# ---------------------------------------------------------------------------

#' Create an anglemaniaObject from a Seurat Object
#'
#' Constructs an \code{\link{anglemaniaObject}} from a given
#' \code{\link[Seurat]{Seurat}} object. This includes extracting and processing
#' count matrices, filtering genes based on expression in a minimum number of
#' cells, and storing results along with dataset and batch information. It also
#' calculates weights for each dataset or batch based on the number of samples.
#'
#' @param seurat_object A \code{\link[Seurat]{Seurat}} object containing
#'   single-cell RNA-seq data.
#' @param dataset_key A character string indicating the column name in the
#'   Seurat object metadata that identifies the dataset to which each cell
#'   belongs. If \code{NA}, all cells are assumed to belong to the same
#'   dataset.
#' @param batch_key A character string indicating the column name(s) in the
#'   Seurat object metadata that identify the batch to which each cell belongs.
#' @param min_cells_per_gene A numeric value indicating the minimum number of
#'   cells in which a gene must be expressed to be included in the analysis.
#'   Default is \code{1}.
#'
#' @return An \code{\link{anglemaniaObject}} containing:
#' \describe{
#'   \item{\code{matrix_list}}{A list of filtered count matrices for each unique
#'     batch.}
#'   \item{\code{dataset_key}}{The dataset key used for splitting the Seurat
#'     object.}
#'   \item{\code{batch_key}}{The batch key used for splitting the Seurat
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
#'   \item Adds a unique batch key to the Seurat object's metadata using
#'     \code{\link{add_unique_batch_key}}.
#'   \item Extracts count matrices for each batch.
#'   \item Filters genes based on the \code{min_cells_per_gene} threshold.
#'   \item Identifies intersected genes present across all batches.
#'   \item Converts count matrices to \code{\link[bigstatsr]{FBM}} objects.
#'   \item Computes weights for each batch or dataset.
#' }
#'
#' @importFrom SeuratObject LayerData
#' @importFrom tidyr unite
#' @importFrom Matrix rowSums
#' @importFrom pbapply pblapply
#' @importFrom dplyr select distinct group_by add_count mutate n_groups
#' @importFrom bigstatsr nb_cores
#' @importFrom checkmate testString
#'
#' @seealso
#' \code{\link{anglemaniaObject-class}},
#' \code{\link{add_unique_batch_key}},
#' \code{\link{anglemania}},
#' \code{\link[bigstatsr]{FBM}}
#' @examples
#' load(system.file(
#'  "extdata",
#'  "seurat_splatter_sim.RData",
#'  package = "anglemania"))
#' batch_key = "Batch"
#' anglemania_object <- create_anglemaniaObject(se,
#'  batch_key = batch_key,
#'  min_cells_per_gene = 1
#'  )
#'  anglemania_object
#' @export create_anglemaniaObject
create_anglemaniaObject <- function(
    seurat_object,
    dataset_key = NA_character_,
    batch_key,
    min_cells_per_gene = 1,
    assay = "RNA",
    allow_missing_features = FALSE,
    min_samples_per_gene = 2
) {
  # Validate inputs
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object needs to be a Seurat object")
  }

  if (checkmate::testString(dataset_key)) {
    if (length(dataset_key) != 1) {
     stop(
       "dataset_key needs to be a character string of length 1 ",
       "corresponding to the column in the metadata of the Seurat ",
       "object that indicates which dataset the cells belong to"
     ) 
    } else if (!(dataset_key %in% colnames(seurat_object[[]]))) {
       stop(
         "dataset_key needs to be a column in the metadata of the Seurat ",
         "object that indicates which dataset the cells belong to"
       )
    }
    message(
      "Using dataset_key: ", dataset_key
    )
  } else if (checkmate::testScalarNA(dataset_key, null.ok = TRUE)) {
      message(
        "No dataset_key specified.\n",
        "Assuming that all samples belong to the same dataset ",
        "and are separated by batch_key: ", batch_key
      )
  }
  else {
    stop(
      "dataset_key needs to be NA/NULL or a character string of length 1 ",
      "corresponding to the column in the metadata of the Seurat ",
      "object that indicates which dataset the cells belong to"
    )
  }

  if (!checkmate::testString(batch_key) || length(batch_key) != 1) {
    stop(
      "batch_key needs to be a character string of length 1 ",
      "corresponding to the column in the metadata of the Seurat ",
      "object that indicates which batch the cells belong to"
    )
  }

  if (!checkmate::testString(assay) || !(assay %in% Assays(seurat_object))) {
    stop(
      "assay needs to be a character string of length 1 ",
      "it needs to correspond to Assays(seurat)"
    )
  }

  # Create unique batch key
  seurat_object <- add_unique_batch_key(
    seurat_object,
    dataset_key,
    batch_key
  )

  meta <- seurat_object[[]]

  # Get the barcodes corresponding to each batch
  # IMPORTANT: When dataset = "string" the order of samples changes
  matrix_list <- split(rownames(meta), meta$batch)
  
  # Extract counts for each batch
  message("Extracting count matrices...")
  message(
    "Filtering each batch to at least ",
    min_cells_per_gene,
    " cells per gene..."
  )
  matrix_list <- lapply(matrix_list, function(cell_barcodes) {
    counts_matrix <- SeuratObject::LayerData(
      seurat_object,
      cells = cell_barcodes,
      layer = "counts",
      assay = assay
    )
    filt_features <- Matrix::rowSums(counts_matrix > 0) >= min_cells_per_gene
    filt_features <- names(filt_features[filt_features])
    counts_matrix <- counts_matrix[filt_features, ]
    return(counts_matrix)
  })

  # checks whether the genes should be present in all samples or just a subset of samples
  if(!allow_missing_features){
    # Reduce to intersection of genes between batches
    message("Using the intersection of filtered genes from all batches...")
    intersect_genes <- Reduce(intersect, lapply(matrix_list, rownames))
  }else{
    message("Using genes which are present in minimally ",min_samples_per_gene, " samples...")
    intersect_genes <- table(unlist(lapply(matrix_list, rownames)))
    intersect_genes <- sort(rownames(intersect_genes)[intersect_genes >= min_samples_per_gene])
  }
  message("Number of genes in intersected set: ", length(intersect_genes))

  
  matrix_list <- pbapply::pblapply(
    matrix_list,
    function(x) {
      genes_in = intersect(rownames(x), intersect_genes)
      x_in <- x[genes_in, ]
      # Checks whether matrix contains missing features
      genes_out = setdiff(intersect_genes, rownames(x))
      if(length(genes_out) > 0){
        # this pads the matrix with missing features
        x_out = Matrix::Matrix(0, nrow=length(genes_out), ncol = ncol(x_in))
        rownames(x_out) = genes_out
        colnames(x_out) = colnames(x_in)
        x_in = rbind(x_in, x_out)
        # IMPORTANT: keeps track that the missing features are ordered in the same way in all matrices
        x_in = x_in[intersect_genes,]
      }
      x_in <- sparse_to_fbm(x_in)
    },
    cl = bigstatsr::nb_cores()
  )

  if (checkmate::testString(dataset_key) && length(unique(meta[[dataset_key]])) > 1) {
    # Calculate the weights only if there are more than two datasets 
    data_info <- meta %>%
      dplyr::select(
        batch,
        all_of(c(dataset_key, batch_key))
      ) %>%
      dplyr::distinct() %>%
      dplyr::group_by(across(all_of(dataset_key))) %>%
      dplyr::add_count(across(all_of(dataset_key)), name = "n_samples") %>%
      dplyr::mutate(
        weight = 1 / n_samples / dplyr::n_groups(.)
      ) %>%
      ## Normalize weights so that their mean is 1 - this helps when the features are not present in all samples
      ungroup() %>%
      mutate(weight = weight / mean(weight))

    weights <- data_info$weight
    names(weights) <- data_info$batch

  } else {
    # if there is only one dataset, weights should be set to one
    data_info <- meta %>%
      dplyr::select(batch, all_of(batch_key)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(weight = 1)

    weights <- data_info$weight
    names(weights) <- data_info$batch
  }
  # Create anglem object
  anglem_object <- new(
    "anglemaniaObject",
    matrix_list = matrix_list,
    dataset_key = ifelse(
      checkmate::testString(dataset_key),
      dataset_key,
      NA_character_
    ),
    batch_key = batch_key,
    data_info = data_info,
    weights = weights,
    min_cells_per_gene = min_cells_per_gene,
    intersect_genes = intersect_genes,
    assay = assay
  )

  return(anglem_object)
}
