#' Identify Empty Droplets and Damaged Cells
#'
#' @name dropletqc_calculation
#'
#' @description This helper function extracts nf meta data and associated cell barcodes from a Seurat object,
#' identifies empty droplets using the `identify_empty_droplets` function, and identifies
#' damaged cells using the `identify_damage_cells` function. The results are added to the
#' Seurat object's meta data.
#'
#' @param Seurat A Seurat object containing the single-cell RNA-seq data.
#'
#' @return The updated Seurat object with droplet QC results added to its meta data.
#'
#' @import Seurat
#'
#' @export
#'
#' @keywords internal

utils::globalVariables(c("nf", "nCount_RNA", "cell_status"))

dropletqc_calculation <- function(Seurat) {

  # Extract nf meta data & associated cell barcode from Seurat object
  edDf <- data.frame(nf = as.numeric(Seurat$nf), umi = Seurat$nCount_RNA)

  # Use droplet_qc function to identify empty droplets
  edresultsDf <- identify_empty_droplets(edDf)

  # Identify damaged cells
  # Create vector of length same as cell number (runs default with no adjustment on groups of cells, need to fill this column requirement while getting results for that sample as a whole)
  n_cells <- length(Cells(Seurat))
  cell_type <- rep(1, n_cells)

  dcDf <- data.frame(
    nf = as.numeric(Seurat$nf),
    umi = Seurat$nCount_RNA,
    cell_status = edresultsDf$cell_status,
    cell_type = cell_type
  )

  # Open a connection to a temporary file for writing
  tmp_conn <- file(tempfile(), open = "wt")

  # Redirect standard output and messages to the temporary file
  sink(tmp_conn)
  sink(tmp_conn, type = "message")

  # Call the (now silenced) verbose function
  dcresultsDf <- identify_damage_cells(dcDf)

  # Reset output redirection
  sink(NULL)
  sink(NULL, type = "message")

  # Close the connection
  close(tmp_conn)

  # Add droplet_qc output to Seurat object meta_data
  Seurat$droplet_qc <- dcresultsDf$df$cell_status

  return(Seurat)
}
