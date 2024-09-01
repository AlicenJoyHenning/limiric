#' limiric
#'
#' @name limiric
#'
#' @description Takes as input filtered scRNA-seq alignment output (or 'Seurat' object)
#' and identifies damaged cells present within the sample according to
#' low dimension mitochondrial and ribosomal clustering.
#'
#' @param project_name String with project or sample name
#' @param filtered_path Directory of filtered alignment output
#' @param seurat_input 'Seurat' object to be used as input over raw files. Default NULL
#' @param min_cells In how many cells should a gene be expressed to be kept
#' @param soupx Perform ambient RNA correction, if TRUE raw_path must be given. Default is FALSE
#' @param raw_path Directory of unfiltered alignment output
#' @param droplet_qc Verify output with droplet_qc, if TRUE velocyto_path must be given. Default is FALSE
#' @param velocyto_path Directory of 'Velocyto' filtered alignment output
#' @param filter_rbc Whether or not red blood cells should be removed. Default is TRUE
#' @param isolate_cd45 Discard non-immune cells. Default is FALSE
#' @param filter_output Should output contain no damaged cells. Default is TRUE
#' @param output_path Directory where 'limiric' output should be generated
#' @param organism "Hsap" if human sample or "Mmus" if mouse sample
#' @param resolution Numeric between 0 and 1.6 describing cluster division. Default 1
#' @param cluster_ranks Numeric describing the number of top ranking clusters to be included as damaged cells. Default 1.
#' @param sample_list Input multiple samples in list. Default is FALSE
#'
#' @return (list) Output storing the final 'Seurat' object
#'
#' @import cowplot
#' @importFrom dplyr %>% mutate case_when
#' @import ggplot2
#' @importFrom methods is
#' @importFrom png readPNG
#' @import Seurat
#' @importFrom utils globalVariables
#'
#' @examples
#'
#' if (interactive()) {
#'
#'   # Load example Seurat object from the limiric package
#'   data("test_data", package = "limiric")
#'
#'
#'   # Run the limiric function with the example data
#'   test <- limiric(
#'     project_name = "test_run",
#'     filter_rbc   = FALSE,
#'     seurat_input = test_data,
#'     output_path  = tempdir()
#'   )
#' }
#'
#' @export

utils::globalVariables(c(
  "current_plots", "droplet_qc", "end_index", "filtered_path", "filter_output",
  "filter_rbc", "i", "isolate_cd45", "min_cells", "nrow", "num_current_plots",
  "num_plots", "organism", "output_path", "page_num", "plots", "plots_per_page",
  "png_path", "project_name", "raw_path", "rds_dir", "rds_files", "results",
  "sample", "sample_list", "seurat_input", "soupx", "temp_result", "test_data",
  "velocyto_path"
))

limiric <- function(

  project_name  = NULL,
  filtered_path = NULL,
  seurat_input  = NULL,
  min_cells     = NULL,
  soupx         = NULL,
  raw_path      = NULL,
  droplet_qc    = NULL,
  velocyto_path = NULL,
  filter_rbc    = NULL,
  isolate_cd45  = NULL,
  filter_output = NULL,
  output_path   = NULL,
  organism      = NULL,
  resolution    = NULL,
  cluster_ranks = NULL,
  sample_list   = NULL

){

  # Two options, either one sample input or multiple samples in the form of a list

  # Option one (default) ------------------------------------

  if (is.null(sample_list)) {

    # Account for defaults
    if (is.null(seurat_input))  {seurat_input = NULL}
    if (is.null(min_cells))     {min_cells = 0}
    if (is.null(resolution))    {resolution = 1}
    if (is.null(cluster_ranks)) {cluster_ranks = 1}
    if (is.null(soupx))         {soupx = FALSE}
    if (is.null(droplet_qc))    {droplet_qc = FALSE}
    if (is.null(filter_rbc))    {filter_rbc = TRUE}
    if (is.null(isolate_cd45))  {isolate_cd45 = FALSE}
    if (is.null(filter_output)) {filter_output = TRUE}
    if (is.null(organism))      {organism = "Hsap"}

    # Run single sample using the limiric_core function
    result <- limiric_core(

      project_name  = project_name,
      filtered_path = filtered_path,
      seurat_input  = seurat_input,
      min_cells     = min_cells,
      resolution    = resolution,
      cluster_ranks = cluster_ranks,
      soupx         = soupx,
      raw_path      = raw_path,
      droplet_qc    = droplet_qc,
      velocyto_path = velocyto_path,
      filter_rbc    = filter_rbc,
      isolate_cd45  = isolate_cd45,
      filter_output = filter_output,
      output_path   = output_path,
      organism      = organism

    )

    # Clean output directories: convert QC rds to png plots

    if (filter_rbc) {
      # Red blood cell QC
      RBCQC_path <- file.path(output_path, "RBCQC", paste0(project_name, "_RBCQC", ".rds"))
      RBCQC      <- readRDS(RBCQC_path)
      ggsave(file.path(output_path, "RBCQC", paste0(project_name, "_RBCQC", ".png")), plot = RBCQC, width = 5, height = 5, dpi = 300)

      file.remove(RBCQC_path)
    }

    if (isolate_cd45) {
      # And immune cell QC
      IMCQC_path <- file.path(output_path, "IMCQC", paste0(project_name, "_IMCQC", ".rds"))
      IMCQC      <- readRDS(IMCQC_path)
      ggsave(file.path(output_path, "IMCQC/", paste0(project_name, "_IMCQC", ".png")), plot = IMCQC, width = 5, height = 5, dpi = 300)

      file.remove(IMCQC_path)
    }

    return(result) # In form of Seurat object

  }


  # Option two ------------------------------------

  else

  {

    # Initialize the output list where resulting Seurat objects will be stored
    results <- list()

    # Run limiric_core function for each sample
    for (i in seq_along(sample_list)) {

      # Define sample
      sample <- sample_list[[i]]

      # Define the inputs from the list
      project_name  <- sample$project_name
      filtered_path <- sample$filtered_path
      seurat_input  <- sample$seurat_input
      min_cells     <- sample$min_cells
      resolution    <- sample$resolution
      cluster_ranks <- sample$cluster_ranks
      soupx         <- sample$soupx
      raw_path      <- sample$raw_path
      droplet_qc    <- sample$droplet_qc
      velocyto_path <- sample$velocyto_path
      filter_rbc    <- sample$filter_rbc
      isolate_cd45  <- sample$isolate_cd45
      filter_output <- sample$filter_output
      output_path   <- sample$output_path
      organism      <- sample$organism

      # Account for defaults
      if (is.null(seurat_input))  {seurat_input = NULL}
      if (is.null(min_cells))     {min_cells = 0}
      if (is.null(resolution))    {resolution = 1}
      if (is.null(cluster_ranks)) {cluster_ranks = 1}
      if (is.null(soupx))         {soupx = FALSE}
      if (is.null(droplet_qc))    {droplet_qc = FALSE}
      if (is.null(filter_rbc))    {filter_rbc = TRUE}
      if (is.null(isolate_cd45))  {isolate_cd45 = FALSE}
      if (is.null(filter_output)) {filter_output = TRUE}
      if (is.null(organism))      {organism = "Hsap"}


      # Call the function with error handling
      # AKA will continue running if any single samples are problematic
      tryCatch({

        temp_result <- limiric_core(

          project_name  = project_name,
          filtered_path = filtered_path,
          seurat_input  = seurat_input,
          min_cells     = min_cells,
          resolution    = resolution,
          cluster_ranks = cluster_ranks,
          soupx         = soupx,
          raw_path      = raw_path,
          droplet_qc    = droplet_qc,
          velocyto_path = velocyto_path,
          filter_rbc    = filter_rbc,
          isolate_cd45  = isolate_cd45,
          filter_output = filter_output,
          output_path   = output_path,
          organism      = organism

        )

        results[[project_name]] <- temp_result

      },

      # If error in one sample, print message and continue with next sample
      error = function(e) {
        cat(paste(e$message, "\n", "Error in processing", project_name, "\n\n"))

      })

    }

    # Condense & clean output directories: convert QC rds to png plots

    if (filter_rbc) {

    # Read in all objects in QC subdirectory of interest
    rds_dir   <- file.path(output_path, "RBCQC")
    rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
    plots     <- lapply(rds_files, readRDS)

    # Loop through the plots & calculate # plot pages needed (max 15 plots each pg)
    num_plots <- length(plots)
    plots_per_page <- 15
    page_num <- 1

    if (num_plots > 0) {
      for (i in seq(1, num_plots, by = plots_per_page)) {
        end_index <- min(i + plots_per_page - 1, num_plots)
        current_plots <- plots[i:end_index]

        # Calculate the number of rows needed for the current set of plots
        num_current_plots <- length(current_plots)
        nrow <- ceiling(num_current_plots / 5)

        # Combine the plots into a single PNG
        png_path <- file.path(output_path, "RBCQC", paste0(project_name, "_RBCQC_", page_num, ".png"))

        create_plot_grid(plots = current_plots,
                         file_path = png_path,
                         nrow = nrow)

        page_num <- page_num + 1

      }

      # Delete all rds files (only want output png)
      file.remove(rds_files)

      }
    }

    if (isolate_cd45) {

      # Read in all objects in QC subdirectory of interest
      rds_dir   <- file.path(output_path, "IMCQC")
      rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
      plots     <- lapply(rds_files, readRDS)

      # Loop through the plots & calculate # plot pages
      num_plots <- length(plots)
      plots_per_page <- 15
      page_num <- 1

      if (num_plots > 0) {
        for (i in seq(1, num_plots, by = plots_per_page)) {
          end_index <- min(i + plots_per_page - 1, num_plots)
          current_plots <- plots[i:end_index]

          # Calculate the number of rows needed for the current set of plots
          num_current_plots <- length(current_plots)
          nrow <- ceiling(num_current_plots / 5)

          # Save the current set of plots as a PNG
          png_path <- file.path(output_path, "IMCQC", paste0(project_name, "_IMCQC_", page_num, ".png"))
          create_plot_grid(plots = current_plots,
                           file_path = png_path,
                           nrow = nrow)

          page_num <- page_num + 1

        }

        # Delete all rds files (only want output png)
        file.remove(rds_files)
      }
    }

    return(results) # list storing each Seurat object

  }

}
