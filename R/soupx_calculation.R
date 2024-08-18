#' Perform SoupX Correction
#'
#' @name soupx_calculation
#'
#' @description This helper function performs SoupX correction with the raw output (must be zipped), estimates the contamination,
#' and outputs the results in a Seurat object for continued workflow.
#'
#' @param raw_path A string representing the path to the raw data.
#' @param table_of_counts A matrix of counts.
#' @param min_cells An integer representing the minimum number of cells.
#' @param project_name A string representing the name of the project.
#'
#' @return A Seurat object with SoupX corrected counts
#'
#' @importFrom cowplot plot_grid
#' @importFrom dplyr %>% mutate
#' @import ggplot2
#' @importFrom png readPNG
#' @import Seurat
#' @import SoupX
#' @importFrom utils globalVariables
#'
#' @export
#'
#' @keywords internal

utils::globalVariables(c("SeuratSoup", "meta_data", "umap_embedding", "adj_matrix", "Seurat"))

soupx_calculation <- function(raw_path,
                              table_of_counts,
                              min_cells,
                              project_name
) {
  # Perform soupx correction with the raw output (must be zipped)
  table_of_droplets <- suppressWarnings(Read10X(raw_path))

  # Run Soup X
  sc <- SoupChannel(table_of_droplets,
                    table_of_counts,
                    calcSoupProfile = FALSE) # Create the soup channel (sc)

  # Estimate the contamination
  sc <- estimateSoup(sc)

  # Use Seurat to cluster the filtered matrix, although not essential it is recommended to get better estimations
  SeuratSoup <- suppressWarnings(CreateSeuratObject(table_of_counts, min.cells = 0))
  SeuratSoup <- suppressWarnings(SCTransform(SeuratSoup, verbose = FALSE) %>%
                                   RunPCA(verbose = FALSE) %>%
                                   RunUMAP(dims = 1:30, verbose = FALSE) %>%
                                   FindNeighbors(dims = 1:30, verbose = FALSE) %>%
                                   FindClusters(verbose = FALSE))

  # Adding the cluster embeddings to the soupx object
  meta_data <- SeuratSoup@meta.data
  umap_embedding <- SeuratSoup@reductions$umap@cell.embeddings
  sc <- suppressWarnings(setClusters(sc, setNames(meta_data$seurat_clusters, rownames(meta_data))))
  sc <- suppressWarnings(setDR(sc, umap_embedding, c("UMAP_1", "UMAP_2")))

  # With defined clusters, run soupx to calculate the contamination fraction rho where rho E (0, 1) and the closer to 1, the more contaminated
  sc <- autoEstCont(sc, verbose = FALSE, doPlot = FALSE)

  # Silencing the output ...
  # Open a connection to a temporary file for writing
  tmp_conn <- file(tempfile(), open = "wt")

  # Redirect standard output and messages to the temporary file
  sink(tmp_conn)
  sink(tmp_conn, type = "message")

  # Call the (now silenced) verbose function
  # Output integer matrix of soup-corrected reads (unzipped output) where contaminated reads are removed
  adj_matrix <- suppressWarnings(adjustCounts(sc, roundToInt = T))

  # Reset output redirection
  sink(NULL)
  sink(NULL, type = "message")

  # Close the connection (end silencing)
  close(tmp_conn)

  # Output results in Seurat object for continued workflow
  Seurat <- suppressWarnings(CreateSeuratObject(counts = adj_matrix, # soupx corrected count matrix
                                                min.cells = min_cells, # At least one cell must express the gene for the gene to be included in the count matrix
                                                min.features = 0,
                                                project = project_name))

  return(Seurat)
}
