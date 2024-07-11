# Class definitions

# -------------------------------------------------------------------------------- #
#' The anglem class
#' 
#' The anglem class is a simple class to construct the correct input for the anglemanise function from a Seurat object and store the results.
#' 
#' @name anglem-class
#' @rdname anglem-class
#' @exportClass anglem
setClass(
    "anglem",
    slots = c(
        matrix_list = "list",
        dataset_key = "character",
        batch_key = "character",
        data_info = "data.frame",
        weights = "numeric",
        list_stats = "list",
        intersect_genes = "character",
        min_cells_per_gene = "numeric",
        integration_genes = "list"),

    prototype = list(
        matrix_list = list(),
        dataset_key = NULL,
        batch_key = NA_character_,
        data_info = data.frame(),
        weights = NA_real_,
        list_stats = list(),
        intersect_genes = NA_character_,
        min_cells_per_gene = 1,
        integration_genes = list(info = "data.frame", genes = "character")
    )
)


# Define a custom show method for anglem class
setMethod("show", "anglem", function(object) {
    cat("Anglem object\n")
    cat("--------------\n")
    cat("Dataset key:", object@dataset_key, "\n")
    cat("Batch key:", object@batch_key, "\n")
    cat("Number of datasets:", ifelse(is.na(object@dataset_key), 1, nrow(unique(object@data_info[, object@dataset_key]))), "\n")
    if(!is.na(object@dataset_key)) cat("Datasets:", paste(unique(object@data_info[, object@dataset_key]), collapse = ", "), "\n")
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
        cat(paste(object@intersect_genes[1:10], collapse = ", "), ", ...\n")
    } else {
        cat(paste(object@intersect_genes, collapse = ", "), "\n")
    }
    cat("Min cells per gene:", object@min_cells_per_gene, "\n")
})


# -------------------------------------------------------------------------------- #
# create an accessor function to get the seurat list and assign a value to it
setGeneric("matrix_list", function(object) standardGeneric("matrix_list"))
setMethod("matrix_list", "anglem", function(object) object@matrix_list)

setGeneric("matrix_list<-", function(object, value) standardGeneric("matrix_list<-"))
setReplaceMethod("matrix_list", "anglem", function(object, value) {
    object@matrix_list <- value
    object
})

# -------------------------------------------------------------------------------- #
# accessor but NO ASSIGNER method for dataset_key
setGeneric("dataset_key", function(object) standardGeneric("dataset_key"))
setMethod("dataset_key", "anglem", function(object) object@dataset_key)

# -------------------------------------------------------------------------------- #
# accessor but NO ASSIGNER method for batch_key
setGeneric("batch_key", function(object) standardGeneric("batch_key"))
setMethod("batch_key", "anglem", function(object) object@batch_key)

# -------------------------------------------------------------------------------- #
# accessor but NO ASSIGNER method for data_info
setGeneric("data_info", function(object) standardGeneric("data_info"))
setMethod("data_info", "anglem", function(object) object@data_info)

# -------------------------------------------------------------------------------- #
# accessor AND assigner method for weights
setGeneric("weights", function(object) standardGeneric("weights"))
setMethod("weights", "anglem", function(object) object@weights)

setGeneric("weights<-", function(object, value) standardGeneric("weights<-"))
setReplaceMethod("weights", "anglem", function(object, value) {
    if (!is.numeric(value)) stop("weights must be numeric")
    if(is.null(names(value))) stop("weight needs to be a named vector")
    # scale the weight so that the sum of weights is 1
    message("scaling provided weights to sum to 1...")
    value = value / sum(value)
    object@weights <- value
    object
})

# -------------------------------------------------------------------------------- #
# accessor and assigner method for list_stats
setGeneric("list_stats", function(object) standardGeneric("list_stats"))
setMethod("list_stats", "anglem", function(object) object@list_stats)

setGeneric("list_stats<-", function(object, value) standardGeneric("list_stats<-"))
setReplaceMethod("list_stats", "anglem", function(object, value) {
    if (!is.list(value)) stop("list_stats must be a list")
    object@list_stats <- value
    object
})

# -------------------------------------------------------------------------------- #
# accessor and assigner method for intersect_genes
setGeneric("intersect_genes", function(object) standardGeneric("intersect_genes"))
setMethod("intersect_genes", "anglem", function(object) object@intersect_genes)

setGeneric("intersect_genes<-", function(object, value) standardGeneric("intersect_genes<-"))
setReplaceMethod("intersect_genes", "anglem", function(object, value) {
    # if (!is.character(value)) stop("intersect_genes must be a character vector")
    object@intersect_genes <- value
    object
})

# -------------------------------------------------------------------------------- #
# accessor and assigner method for integration_genes
setGeneric("extract_integration_genes", function(object) standardGeneric("extract_integration_genes"))
setMethod("extract_integration_genes", "anglem", function(object){
    object@integration_genes$genes
})


# -------------------------------------------------------------------------------- #
#' Add a unique batch key to the metadata of a seurat object
#' 
#' @param seurat_object A Seurat object
#' @param dataset_key A column name in the metadata of the Seurat object
#' @param batch_key A column name in the metadata of the Seurat object
#' @return A Seurat object with a unique batch key
#' @export add_unique_batch_key
#' @examples
#' add_unique_batch_key(seurat_object, dataset_key = "dataset", batch_key = "batch")
add_unique_batch_key <- function(seurat_object, dataset_key = NULL, batch_key, new_unique_batch_key = "batch") {

    if (!is.null(dataset_key)) {
                meta <- seurat_object[[]] %>%
                    tidyr::unite("batch", all_of(batch_key), sep = "_", remove = FALSE) %>%
                    tidyr::unite("batch", all_of(c(dataset_key, "batch")), sep = ":", remove = FALSE)
            } else {
                meta <- seurat_object[[]] %>%
                    tidyr::unite("batch", all_of(batch_key), sep = "_", remove = FALSE)
            }
    seurat_object[[]] <- meta

    return(seurat_object)
}

# -------------------------------------------------------------------------------- #

#' Create an Anglem Object from a Seurat Object
#' 
#' This function constructs an `anglem` object from a given Seurat object, which includes the extraction and processing of count matrices, filtering of genes based on expression in a minimum number of cells, and storing the results along with dataset and batch information.
#' It also adds weights for each dataset based on the number of batches in each dataset. If n(dataset) == 1, then the weight is 1/n(batches) for each batch.
#' 
#' @param seurat_object A Seurat object containing single-cell RNA-seq data.
#' @param dataset_key A character vector of length 1 indicating the column name in the Seurat object meta data that identifies the dataset to which each cell belongs. If NULL, all cells are assumed to belong to the same dataset.
#' @param batch_key A character vector indicating the column name(s) in the Seurat object meta data that identify the batch to which each cell belongs.
#' @param min_cells_per_gene A numeric value indicating the minimum number of cells in which a gene must be expressed to be included in the analysis. Default is 10.
#' 
#' @return An `anglem` object containing:
#' \item{matrix_list}{A list of filtered count matrices, one for each unique combination of dataset and batch.}
#' \item{dataset_key}{The dataset key used for splitting the Seurat object.}
#' \item{batch_key}{The batch key used for splitting the Seurat object.}
#' \item{data_info}{A data frame summarizing the number of samples per dataset and their weights.}
#' \item{weights}{A numeric vector of weights for each dataset based on the number of samples.}
#' \item{intersect_genes}{A character vector of genes that are expressed in at least the specified number of cells across all batches.}
#' \item{min_cells_per_gene}{The minimum number of cells per gene threshold used for filtering.}
#' 
#' @importFrom Seurat SplitObject
#' @importFrom Seurat LayerData
#' @importFrom tidyr unite
#' @importFrom Matrix Matrix
#' @import dplyr
#' 
#' @export create_anglem
#' 
#' @examples
#' \dontrun{
#' se <- CreateSeuratObject(counts = matrix(rpois(2000, lambda = 10), nrow = 200, ncol = 10))
#' se$Experiment <- rep(c("Exp1", "Exp2"), each = 5)
#' se$Method <- rep(c("M1", "M2"), times = 5)
#' test <- create_anglem(se, dataset_key = "Experiment", batch_key = "Method")
#' }
# constructor for anglem class
create_anglem <- function(
    seurat_object,
    dataset_key = NULL,
    batch_key,
    min_cells_per_gene = 1
    ){
        # validate inputs
        if (class(seurat_object) != "Seurat") {
            stop("seurat_object needs to be a Seurat object")
        }

        if (is.null(dataset_key)) {
            message("No dataset_key specified.\n Assuming that all samples belong to the same dataset and are separated by batch_key: ", batch_key)
        }

        if (!is.null(dataset_key)) {
            if (!is.character(dataset_key) || length(dataset_key) != 1) {
                stop("dataset_key needs to be a character vector of length 1 corresponding to the column in the meta data of the Seurat object that indicates which dataset the cells belong to")
            }
        }

        if (is.null(batch_key)){
            stop("batch_key is not specified! batch_key needs to be a character vector indicating the batch key(s) in the meta data of the Seurat object")
        }

        if (!is.character(batch_key)) {
            stop("batch_key needs to be a character vector indicating the batch key(s) in the meta data of the Seurat object")
        }

        # create column in meta data that combines dataset_key and batch_key for each cell so that you can split the Seurat Object by this column
        seurat_object <- add_unique_batch_key(seurat_object, dataset_key, batch_key)

        meta <- seurat_object[[]]
        message("Extracting count matrices...") 
        matrix_list <- Seurat::SplitObject(seurat_object, split.by = "batch")
        names(matrix_list) <- unique(meta$batch)
        matrix_list <- lapply(matrix_list, function(x) {
            x <- SeuratObject::LayerData(x, layer = "counts", assay = "RNA") # GetAssayData or LayerData from SeuratObject?
                })

        if (!is.null(dataset_key)) {
            data_info <- meta %>%
                dplyr::select(batch, all_of(c(dataset_key, batch_key))) %>%
                dplyr::distinct() %>%
                dplyr::group_by(across(all_of(dataset_key))) %>%
                dplyr::add_count(across(all_of(dataset_key)), name = "n_samples") %>%
                dplyr::mutate(weight = 1/n_samples/n_groups(.)) # equal weight for each sample of a dataset. ==> sum(weight) = 1

            weights <- data_info$weight
            names(weights) <- data_info$batch
        } else { # if there is only one dataset, then weight is 1 ==> no need to adjust. 
            data_info <- meta %>%
                dplyr::select(batch, all_of(batch_key)) %>%
                dplyr::distinct() %>%
                dplyr::mutate(weight = 1/nrow(.)) 

            weights <- data_info$weight
            names(weights) <- data_info$batch
        }
        # create anglem object
        anglem_object <- new(
            "anglem",
            matrix_list = matrix_list, # matrix_list: A list of sparse matrices, one for each batch.
            dataset_key = ifelse(is.null(dataset_key), NA_character_, dataset_key), # dataset_key: The key used to denote the dataset.
            batch_key = batch_key, # batch_key: The key used to denote the batch.
            data_info = data_info, # data_info: A data frame summarizing the number of samples per dataset and their weights.
            weights = weights, # weights: A numeric vector of weights for each dataset based on the number of samples. By default, the weights are adjusted to the size of the dataset_key if provided. Otherwise each batch has the same weight.
            min_cells_per_gene = min_cells_per_gene # min_cells_per_gene: The minimum number of cells per gene threshold used for filtering.
        )

        message("Filtering each batch to at least ", min_cells_per_gene, " cells per gene...")
        intersect_genes(anglem_object) <- Reduce(intersect, lapply(matrix_list(anglem_object), function(x){
            # indicate which features are expressed 
            # then only keep the genes that are expressed in at least min_cells_per_gene cells
            filt_features <- Matrix::rowSums(x > 0) >= anglem_object@min_cells_per_gene  

            names(filt_features[filt_features])
        }))

        message("Using the intersection of filtered genes from all batches...")
        message("Number of genes in intersected set: ", length(intersect_genes(anglem_object)))

        matrix_list(anglem_object) <- pbapply::pblapply(matrix_list(anglem_object), function(x){
            x <- x[intersect_genes(anglem_object), ]
            x <- sparse_to_fbm(x) # convert sparse matrix into FBM. This function is part of anglemania and is described in anglemanise-utils.R
        }, cl = 4) # COMMENT: add cores parameter??

        return(anglem_object)
}

