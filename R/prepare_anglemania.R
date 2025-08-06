# ---------------------------------------------------------------------------
#' @describeIn anglemania Check Parameters provided to the anglemania function
#' @import checkmate
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @examples
#' sce <- sce_example()
#' params <- check_params(
#'   sce,
#'   batch_key = "batch",
#'   dataset_key = "dataset",
#'   zscore_mean_threshold = 2.5,
#'   zscore_sn_threshold = 2.5,
#'   max_n_genes = 2000,
#'   method = "cosine",
#'   min_cells_per_gene = 1,
#'   min_samples_per_gene = 2,
#'   allow_missing_features = FALSE,
#'   permute_row_or_column = "column",
#'   permutation_function = "sample",
#'   prefilter_threshold = 0.5,
#'   normalization_method = "divide_by_total_counts",
#'   verbose = TRUE
#' )
#' @return A list of validated parameters
#' @export
check_params <- function(
    sce,
    batch_key,
    dataset_key,
    zscore_mean_threshold,
    zscore_sn_threshold,
    max_n_genes,
    method,
    min_cells_per_gene,
    min_samples_per_gene,
    allow_missing_features,
    permute_row_or_column,
    permutation_function,
    prefilter_threshold,
    normalization_method,
    verbose
) {
    # Validate inputs
    # sce
    checkmate::assert_class(
        sce,
        c("SingleCellExperiment", "SummarizedExperiment")
    )
    # batch_key
    if (
        !checkmate::test_string(batch_key) ||
            length(batch_key) != 1 ||
            !(batch_key %in% colnames(colData(sce)))
    ) {
        stop(
            "batch_key needs to be a character string of length 1 ",
            "corresponding to the column in the metadata of the SingleCellExperiment ",
            "object that indicates which batch the cells belong to"
        )
    }
    # dataset_key
    if (checkmate::test_string(dataset_key)) {
        if (length(dataset_key) != 1) {
            stop(
                "dataset_key needs to be a character string of length 1 ",
                "corresponding to the column in the metadata of the Seurat ",
                "object that indicates which dataset the cells belong to"
            )
        } else if (!(dataset_key %in% colnames(colData(sce)))) {
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
        # if no dataset key is provided, make it NA, because the anglemania
        # object class expects it to be NA or a character
        dataset_key <- NA_character_
    } else {
        stop(
            "dataset_key needs to be NA/NULL or a character string of length 1 ",
            "corresponding to the column in the metadata of the Seurat ",
            "object that indicates which dataset the cells belong to"
        )
    }
    # zscore_mean_threshold
    checkmate::assert_numeric(zscore_mean_threshold, lower = 0, len = 1)
    # zscore_sn_threshold
    checkmate::assert_numeric(zscore_sn_threshold, lower = 0, len = 1)
    # max_n_genes
    checkmate::assert_integerish(
        max_n_genes,
        lower = 1,
        len = 1,
        null.ok = TRUE
    )
    # method
    checkmate::assert_choice(
        method,
        c("cosine", "spearman")
    )
    # min_cells_per_gene
    checkmate::assert_integerish(min_cells_per_gene, lower = 1, len = 1)
    # min_samples_per_gene
    checkmate::assert_integerish(min_samples_per_gene, lower = 1, len = 1)
    # allow_missing_features
    checkmate::assert_logical(allow_missing_features, len = 1)
    # permute_row_or_column
    checkmate::assert_choice(
        permute_row_or_column,
        c("row", "column")
    )
    # permutation_function
    checkmate::assert_choice(
        permutation_function,
        c("sample", "permute_nonzero")
    )
    # prefilter_threshold
    checkmate::assert_numeric(prefilter_threshold, lower = 0, len = 1)
    # normalization_method
    checkmate::assert_choice(
        normalization_method,
        c("divide_by_total_counts", "scale_by_total_counts")
    )
    # verbose
    checkmate::assert_logical(
        verbose,
        len = 1
    )
    return(list(
        batch_key = batch_key,
        dataset_key = dataset_key,
        zscore_mean_threshold = zscore_mean_threshold,
        zscore_sn_threshold = zscore_sn_threshold,
        max_n_genes = max_n_genes,
        method = method,
        min_cells_per_gene = min_cells_per_gene,
        min_samples_per_gene = min_samples_per_gene,
        allow_missing_features = allow_missing_features,
        permute_row_or_column = permute_row_or_column,
        permutation_function = permutation_function,
        prefilter_threshold = prefilter_threshold,
        normalization_method = normalization_method,
        verbose = verbose
    ))
}


#' Internal Helper function to calculate and set weights for the batches
#' in a SCE object
#' @param col_data A data frame containing the metadata of the SCE/SE object
#' @param batch_key A character string specifying the column name in the metadata that identifies the batch
#' @param dataset_key A character string specifying the column name in the metadata that identifies the dataset
#' @import dplyr
#' @importFrom bigstatsr nb_cores
#' @importFrom checkmate test_string
#' @import dplyr
#' @keywords internal
#' @noRd
.set_weights <- function(
    col_data,
    batch_key,
    dataset_key = NA_character_
) {
    col_data <- as.data.frame(col_data)
    if (
        checkmate::test_string(dataset_key) &&
            length(unique(col_data[[dataset_key]])) > 1
    ) {
        data_info <- col_data |>
            dplyr::select(
                anglemania_batch,
                dplyr::all_of(c(dataset_key, batch_key))
            ) |>
            dplyr::distinct() |>
            dplyr::group_by(dplyr::across(dplyr::all_of(dataset_key))) |>
            dplyr::add_count(name = "n_samples") |>
            dplyr::mutate(weight = 1 / n_samples / dplyr::n_groups(.)) |>
            dplyr::ungroup() |>
            dplyr::mutate(weight = weight / mean(weight))
    } else {
        data_info <- col_data |>
            dplyr::select(anglemania_batch, dplyr::all_of(batch_key)) |>
            dplyr::distinct() |>
            dplyr::mutate(weight = 1)
    }
    data_info
}


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
#' @import dplyr
#' @import checkmate
#' @examples
#' sce <- sce_example()
#' head(SummarizedExperiment::colData(sce))
#' sce <- add_unique_batch_key(
#'     sce,
#'     batch_key = "batch",
#'     dataset_key = "dataset"
#' )
#' head(SummarizedExperiment::colData(sce))
#' @describeIn anglemania Temporarily add a unique batch key to the dataset
#' @export
add_unique_batch_key <- function(
    sce,
    dataset_key = NA_character_,
    batch_key
) {
    col_data <- SummarizedExperiment::colData(sce) |> as.data.frame()
    if (checkmate::test_string(dataset_key)) {
        col_data <- col_data |>
            tidyr::unite(
                "batch",
                dplyr::all_of(batch_key),
                sep = "_",
                remove = FALSE
            ) |>
            tidyr::unite(
                "batch",
                dplyr::all_of(c(dataset_key, "batch")),
                sep = ":",
                remove = FALSE
            )
    } else {
        col_data <- col_data |>
            tidyr::unite(
                "batch",
                dplyr::all_of(batch_key),
                sep = "_",
                remove = FALSE
            )
    }

    SummarizedExperiment::colData(sce)$anglemania_batch <- col_data$batch
    return(sce)
}


#' Internal Helper function to print messages
#' @param verbose Logical indicating whether to print messages
#' @param ... Message to print
#' @keywords internal
#' @noRd
vmessage <- function(verbose, ...) {
    if (isTRUE(verbose)) message(...)
}


prepare_matrices <- function(
    matrix_list,
    intersect_genes,
    verbose = TRUE
) {
    #------ PADDING MATRICES & FBM CONVERSION ------#
    # Extract counts for each batch
    vmessage(verbose, "Extracting count matrices...")

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
                # IMPORTANT: keeps track that the missing features are
                # ordered in the same way in all matrices
            } else {
                x <- x[match(intersect_genes, rownames(x)), ]
            }
            sparse_to_fbm(x)
        },
        cl = min(4, bigstatsr::nb_cores())
    )
    return(matrix_list)
}

#'
#' @describeIn anglemania Extract the intersected genes from a
#' list of matrices (count matrices from different batches/datasets).
#' It also allows for missing features in individual matrices, so
#' that a feature does not have to be present in every single batch.
#' @param matrix_list A list of \code{bigstatsr::FBM} objects.
#' @param allow_missing_features Logical indicating whether
#' to allow missing features.
#' @param min_samples_per_gene Integer indicating the minimum
#' number of samples per gene.
#' @param verbose Logical indicating whether to print messages.
#' @examples
#' library(SingleCellExperiment)
#' sce <- sce_example()
#' barcodes_by_batch <- split(colnames(sce), colData(sce)$batch)
#' matrix_list <- lapply(barcodes_by_batch, function(barcodes) {
#'     SingleCellExperiment::counts(sce)[, barcodes]
#' })
#' intersect_genes <- get_intersect_genes(matrix_list)
#' head(intersect_genes)
#' @return A character vector of intersected genes.
#' @export
get_intersect_genes <- function(
    matrix_list,
    allow_missing_features = FALSE,
    min_samples_per_gene = 1,
    verbose = TRUE
) {
    #------ INTERSECT GENES ------#
    # checks whether the genes should be present in all samples or
    # just a subset of samples
    if (!allow_missing_features) {
        # Reduce to intersection of genes between batches
        vmessage(
            verbose,
            "Using the intersection of filtered genes from all batches..."
        )
        intersect_genes <- Reduce(intersect, lapply(matrix_list, rownames))
    } else {
        vmessage(
            verbose,
            "Using genes which are present in minimally ",
            min_samples_per_gene,
            " samples..."
        )
        intersect_genes <- table(unlist(lapply(matrix_list, rownames)))
        intersect_genes <- sort(rownames(intersect_genes)[
            intersect_genes >= min_samples_per_gene
        ])
    }
    vmessage(
        verbose,
        "Number of genes in intersected set: ",
        length(intersect_genes)
    )
    return(intersect_genes)
}
