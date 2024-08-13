#' limiric
#'
#' Function that takes as input filtered scRNA-seq alignment output (or a **Seurat** object)
#' and identifies damaged cells present within the sample according to
#' low dimension mitochondrial and ribosomal clustering.
#'
#' @name limiric
#'
#' @param ProjectName String with project or sample name
#' @param FilteredPath Directory of filtered alignment output
#' @param SeuratInput Seurat object to be used as input over raw files. Default NULL
#' @param MinCells In how many cells should a gene be expressed to be kept
#' @param SoupX Perform ambient RNA correction, if TRUE RawPath must be given. Default is FALSE
#' @param RawPath Directory of unfiltered alignment output
#' @param DropletQC Verify output with DropletQC, if TRUE VelocytoPath must be given. Default is FALSE
#' @param VelocytoPath Directory of Veocyto filtered alignment output
#' @param FilterRBC Whether or not red blood cells should be removed. Default is TRUE
#' @param IsolateCD45 Discard non-immune cells. Default is FALSE
#' @param FilterOutput Should output contain no damaged cells. Default is TRUE
#' @param OutputPath Directory where limiric output should be generated
#' @param Organism "Hsap" if human sample or "Mmus" if mouse sample
#' @param sample_list Input multiple samples in list. Default is FALSE
#'
#' @return (list) Output storing the final **Seurat** object
#'
#' @import cowplot
#' @importFrom dplyr %>% pull group_by summarise mutate arrange slice case_when
#' @import ggplot2
#' @import magick
#' @import Matrix
#' @import png
#' @import Seurat
#' @import SoupX
#' @importFrom stats setNames
#' @importFrom utils data write.csv
#'
#'
#' @examples
#' \dontrun{
#' ## 1. Basic usage
#'
#' # Detect damaged cells in a single sample using
#' # the filtered alignment files
#'
#' SRR1234567 <- limiric(
#'     ProjectName  = "SRR1234567",
#'     FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
#'     OutputPath   = "/home/user/alignment/limiric/"
#' )
#'
#'
#' ## 2. Perform ambient RNA correction prior to damaged cell detection
#'
#' # Detect damaged cells after performing ambient RNA correction
#'
#' SRR1234567 <- limiric(
#'     ProjectName  = "SRR1234567",
#'     FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
#'     SoupX        = TRUE,
#'     RawPath      = "/home/user/alignment/SRR1234567/raw/",
#'     OutputPath   = "/home/user/alignment/limiric/"
#' )
#'
#'
#' ## 3. First isolate the immune cells present in the sample
#' # then identify damaged cells
#'
#' SRR1234567 <- limiric(
#'     ProjectName  = "SRR1234567",
#'     FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
#'     IsolateCD45  = TRUE,
#'     OutputPath   = "/home/user/alignment/limiric/"
#' )
#'
#'
#' ## 4. Combine limiric damaged cell annotations with DropletQC
#'
#' # Detect damaged cells and compare results with those from DropletQC
#'
#' SRR1234567 <- limiric(
#'     ProjectName  = "SRR1234567",
#'     FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
#'     DropletQC    = TRUE,
#'     VelocytoPath = "/home/user/alignment/velocyto/",
#'     OutputPath   = "/home/user/alignment/limiric/"
#' )
#'
#'
#' ## 5. Combine the previous four conditions
#'
#' SRR1234567 <- limiric(
#'     ProjectName  = "SRR1234567",
#'     FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
#'     SoupX        = TRUE,
#'     RawPath      = "/home/user/alignment/SRR1234567/raw/",
#'     DropletQC    = TRUE,
#'     IsolateCD45  = TRUE,
#'     VelocytoPath = "/home/user/alignment/velocyto/",
#'     OutputPath   = "/home/user/alignment/limiric/"
#' )
#'
#'
#' ## 6. Conditions exactly as for 5
#' # but multiple samples are being run instead of just one
#'
#' sample_list <- list(
#'
#'     list(ProjectName  = "SRR1234567",
#'          FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
#'          SoupX        = TRUE,
#'          RawPath      = "/home/user/alignment/SRR1234567/raw/",
#'          DropletQC    = TRUE,
#'          IsolateCD45  = TRUE,
#'          VelocytoPath = "/home/user/alignment/velocyto/",
#'          OutputPath   = "/home/user/alignment/limiric/"),
#'
#'     list(ProjectName  = "SRR1234568",
#'          FilteredPath = "/home/user/alignment/SRR1234568/filtered/",
#'          SoupX        = TRUE,
#'          RawPath      = "/home/user/alignment/SRR1234568/raw/",
#'          DropletQC    = TRUE,
#'          IsolateCD45  = TRUE,
#'          VelocytoPath = "/home/user/alignment/velocyto/",
#'          OutputPath   = "/home/user/alignment/limiric/"),
#'
#'     list(ProjectName  = "SRR1234569",
#'          FilteredPath = "/home/user/alignment/SRR1234569/filtered/",
#'          SoupX        = TRUE,
#'          RawPath      = "/home/user/alignment/SRR1234569/raw/",
#'          DropletQC    = TRUE,
#'          IsolateCD45  = TRUE,
#'          VelocytoPath = "/home/user/alignment/velocyto/",
#'          OutputPath   = "/home/user/alignment/limiric/")
#' )
#'
#' GSE1234567 <- limiric(sample_list = sample_list)
#'
#' ## 7. Alternatively, use a Seurat object as input
#'
#' SRR1234567 <- limiric(
#'     ProjectName  = "SRR1234567",
#'     SeuratInput  = seurat_object,
#'     OutputPath   = "/home/user/alignment/limiric/"
#' )
#' }
#'
#' @export

utils::globalVariables(c(
  "OutputPath", "ProjectName", "rds_dir", "rds_files", "plots", "num_plots",
  "plots_per_page", "page_num", "i", "end_index", "current_plots", "num_current_plots",
  "nrow", "png_path", "sample_list", "results", "sample", "FilteredPath",
  "SeuratInput", "MinCells", "SoupX", "RawPath", "DropletQC", "VelocytoPath",
  "FilterRBC", "IsolateCD45", "FilterOutput", "Organism", "temp_result"
))

limiric <- function(

  ProjectName  = NULL,     # String with the name of the project
  FilteredPath = NULL,     # Path to zipped filtered output of alignment
  SeuratInput  = NULL,     # Whether or not the input is a seurat object
  MinCells     = NULL,     # Minimum number of cells a gene must be expressed in to be retained
  SoupX        = NULL,     # Default is not to do SoupX ambient RNA filtering
  RawPath      = NULL,     # Path to zipped raw output of alignment
  DropletQC    = NULL,     # Default is not to include nuclear fraction data
  VelocytoPath = NULL,     # Path to zipped filtered output of velocyto
  FilterRBC    = NULL,     # Default is to remove red bloods cells
  IsolateCD45  = NULL,     # default is not to filter immune cells
  FilterOutput = NULL,     # Default is to remove damaged cells
  OutputPath   = NULL,     # default output directory is the working directory
  Organism     = NULL,
  sample_list  = NULL

){

  # First, handle two possible cases, either one sample as input or multiple in the form of a list

  if (is.null(sample_list)) {

    # Account for defaults
    if (is.null(SeuratInput)) { SeuratInput = NULL}
    if (is.null(MinCells)) { MinCells = 0 }
    if (is.null(SoupX)) { SoupX = FALSE}
    if (is.null(DropletQC)) { DropletQC = FALSE}
    if (is.null(FilterRBC)) { FilterRBC = TRUE}
    if (is.null(IsolateCD45)) { IsolateCD45 = FALSE}
    if (is.null(FilterOutput)) { FilterOutput = TRUE}
    if (is.null(Organism)) { Organism = "Hsap"}

    # Run single sample (default) using the core function
    result <- limiricCore(

      ProjectName  = ProjectName,
      FilteredPath = FilteredPath,
      SeuratInput  = SeuratInput,
      MinCells     = MinCells,
      SoupX        = SoupX,
      RawPath      = RawPath,
      DropletQC    = DropletQC,
      VelocytoPath = VelocytoPath,
      FilterRBC    = FilterRBC,
      IsolateCD45  = IsolateCD45,
      FilterOutput = FilterOutput,
      OutputPath   = OutputPath,
      Organism     = Organism

    )


    # Clean output directories: convert QC rds to png plots

    if (FilterRBC) {

    # Red blood cell QC
    RBCQC_path <- file.path(OutputPath, "RBCQC", paste0(ProjectName, "_RBCQC", ".rds"))
    RBCQC      <- readRDS(RBCQC_path)
    ggsave(file.path(OutputPath, "RBCQC", paste0(ProjectName, "_RBCQC", ".png")), plot = RBCQC, width = 5, height = 5, dpi = 300)

    file.remove(RBCQC_path)

    }

    if (IsolateCD45 == TRUE) {

      # And immune QC
      IMCQC_path <- file.path(OutputPath, "IMCQC", paste0(ProjectName, "_IMCQC", ".rds"))
      IMCQC      <- readRDS(IMCQC_path)
      ggsave(file.path(OutputPath, "IMCQC/", paste0(ProjectName, "_IMCQC", ".png")), plot = IMCQC, width = 5, height = 5, dpi = 300)
      file.remove(IMCQC_path)

    }

    return(result)

  }

  # If there is a sample list, extract elements from each list item and run with core PreProcess function
  else

  {

    # Initialize the output list where resulting Seurat objects will be stored
    results <- list()

    # Run PreProcess function for each sample
    for (i in seq_along(sample_list)) {

      # Define sample
      sample <- sample_list[[i]]

      # Define the inputs from the list
      ProjectName  <- sample$ProjectName
      FilteredPath <- sample$FilteredPath
      SeuratInput  <- sample$SeuratInput
      MinCells     <- sample$MinCells
      SoupX        <- sample$SoupX
      RawPath      <- sample$RawPath
      DropletQC    <- sample$DropletQC
      VelocytoPath <- sample$VelocytoPath
      FilterRBC    <- sample$FilterRBC
      IsolateCD45  <- sample$IsolateCD45
      FilterOutput <- sample$FilterOutput
      OutputPath   <- sample$OutputPath
      Organism     <- sample$Organism

      # Account for defaults
      if (is.null(SeuratInput)) { SeuratInput = NULL}
      if (is.null(MinCells)) { MinCells = 0 }
      if (is.null(SoupX)) { SoupX = FALSE}
      if (is.null(DropletQC)) { DropletQC = FALSE}
      if (is.null(FilterRBC)) { FilterRBC = TRUE}
      if (is.null(IsolateCD45)) { IsolateCD45 = FALSE}
      if (is.null(FilterOutput)) { FilterOutput = TRUE}
      if (is.null(Organism)) { Organism = "Hsap"}

      # Call the core function with error handling
      tryCatch({

        temp_result <- limiricCore(

          ProjectName  = ProjectName,
          FilteredPath = FilteredPath,
          SeuratInput  = SeuratInput,
          MinCells     = MinCells,
          SoupX        = SoupX,
          RawPath      = RawPath,
          DropletQC    = DropletQC,
          VelocytoPath = VelocytoPath,
          FilterRBC    = FilterRBC,
          IsolateCD45  = IsolateCD45,
          FilterOutput = FilterOutput,
          OutputPath   = OutputPath,
          Organism     = Organism

        )

        results[[ProjectName]] <- temp_result

      },

      # If error in one sample, print message and continue with next sample
      error = function(e) {
        cat(paste(e$message, "\n", "Error in processing", ProjectName, "\n", "Likely low cell number\n\n"))

      })

    }


    if (FilterRBC == TRUE) {

    # Collate output QC plots ####
    # RBC
    rds_dir   <- file.path(OutputPath, "RBCQC") # Define the directory containing the RDS files
    rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
    plots     <- lapply(rds_files, readRDS) # Read all RDS files into a list of ggplot objects

    # Function to save a grid of plots(5 in each row and max 3 rows)
    save_plot_grid <- function(plots, file_path, ncol = 5, nrow = 3) {
      plot_grid <- plot_grid(plotlist = plots, ncol = ncol, nrow = nrow)
      ggsave(file_path, plot = plot_grid, width = ncol * 5, height = nrow * 5, dpi = 300)
    }

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

        # Save the current set of plots as a PNG
        png_path <- file.path(OutputPath, "RBCQC", paste0(ProjectName, "_RBCQC_", page_num, ".png"))
        save_plot_grid(current_plots, png_path, ncol = 5, nrow = nrow)

        page_num <- page_num + 1
      }

      # Delete all rds files (only want output image)
      file.remove(rds_files)
    }

    }

    if (!is.null(IsolateCD45) & IsolateCD45 == TRUE) {

      # Collate output IMC QC plots ####
      # IMC
      rds_dir   <- file.path(OutputPath, "IMCQC") # Define the directory containing the RDS files
      rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
      plots     <- lapply(rds_files, readRDS) # Read all RDS files into a list of ggplot objects

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
          png_path <- file.path(OutputPath, "IMCQC", paste0(ProjectName, "_IMCQC_", page_num, ".png"))
          save_plot_grid(current_plots, png_path, ncol = 5, nrow = nrow)

          page_num <- page_num + 1
        }

        # Delete all rds files (only want output image)
        file.remove(rds_files)
      }
    }

    return(results)

  }

}
