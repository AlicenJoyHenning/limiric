#' test_data
#'
#' Small Seurat object used for demonstrating the `limiric` function.
#'
#' @format A Seurat object with the following structure:
#' \describe{
#'   \item{assays}{List of assays, each containing a matrix of gene expression data.}
#'   \item{meta.data}{Data frame containing metadata for each cell.}
#'   \item{reductions}{List of dimensionality reduction results.}
#' }
#' @source Simulated data
#' @examples
#' data(test_data)
#' # Detect damaged cells in the Seurat object
#' result <- limiric(
#'   project_name  = "ExampleProject",
#'   seurat_input  = test_data,
#'   output_path   = tempdir()
#' )
#' # Check the result
#' print(result)
"test_data"
