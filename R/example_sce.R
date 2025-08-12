#' Generate example SingleCellExperiment object
#'
#' This function generates a SingleCellExperiment object with 2 batches and 2
#' datasets. The object contains 300 genes and 600 cells. The counts matrix is
#' generated using the \code{rpois} function with different lambda values
#' for the two batches.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats rpois
#' @importFrom S4Vectors DataFrame
#' @importFrom withr with_seed
#' @param seed The seed for the random number generator. Because this function
#' is only used locally, we allow to set a seed within the function.
#' @return A SingleCellExperiment object
#' @examples
#' sce <- sce_example()
#' sce
#' @export
sce_example <- function(seed = 42) {
    # generate example data with two "batches"
    withr::with_seed(seed, {
        counts_mat <- cbind(
            matrix(
                rpois(300 * 300, lambda = 5),
                nrow = 300,
                dimnames = list(
                    paste0("gene", seq_len(300)),
                    paste0("cell", seq_len(300))
                )
            ),
            matrix(
                rpois(300 * 300, lambda = 3),
                nrow = 300,
                dimnames = list(
                    paste0("gene", seq_len(300)),
                    paste0("cell", seq_len(300) + 300)
                )
            )
        )
    })
    # add metadata
    batch <- rep(c("batch1", "batch2"), each = 300)
    dataset <- rep(
        c("dataset1", "dataset2", "dataset1", "dataset2"),
        each = 150
    )
    col_data <- S4Vectors::DataFrame(batch = batch, dataset = dataset)

    # Create the SingleCellExperiment object
    sce <- SingleCellExperiment(
        assays = list(counts = counts_mat),
        colData = col_data
    )
    return(sce)
}
