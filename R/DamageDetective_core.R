#' DamageDetective_core
#'
#' @name DamageDetective_core
#'
#' @description Core helper function enabling the DamageDetective modules to run
#'
#' @param project_name String with project or sample name
#' @param filtered_path Directory of filtered alignment output
#' @param seurat_input Seurat object to be used as input over raw files. Default NULL
#' @param min_cells In how many cells should a gene be expressed to be kept
#' @param soupx Perform ambient RNA correction, if TRUE raw_path must be given. Default is FALSE
#' @param raw_path Directory of unfiltered alignment output
#' @param filter_rbc Whether or not red blood cells should be removed. Default is TRUE
#' @param hemo_threshold Percent hemoglobin expression above which cells are filtered. Default is 20
#' @param isolate_cd45 Discard non-immune cells. Default is FALSE
#' @param filter_output Should output contain no damaged cells. Default is TRUE
#' @param output_path Directory where output can be generated
#' @param resolution Numeric between 0 and 1.6 describing cluster division. Default 1
#' @param cluster_ranks Numeric describing the number of top ranking clusters to be included as damaged cells. Default 1.
#' @param organism "Hsap" if human sample or "Mmus" if mouse sample
#' @param verbose Print messages to the console. Default is TRUE
#'
#' @return (list) Output list storing the final Seurat object
#'
#' @import cowplot
#' @importFrom dplyr %>% mutate case_when
#' @import ggplot2
#' @importFrom png readPNG
#' @import Seurat
#' @importFrom utils write.csv
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Load test data
#' data("test_data", package = "DamageDetective")
#'
#' test <- DamageDetective(
#'   project_name  = "test",
#'   seurat_input  = test_data,
#'   filter_rbc    = FALSE,
#'   output_path   = tempdir()
#' )
#'
#' }
#'
#' @keywords internal

utils::globalVariables(c(
  "avg_mt_percent", "avg_rb_percent", "combined_rank", "gene_name",
  "hemo.percent", "human_annotations", "MtRb", "mouse_annotations",
  "mt.percent", "nCount_RNA", "nFeature_RNA", "ptprc.percent", "QC",
  "quality", "rank_mt", "rank_rb", "rb.percent", "seurat_clusters",
  "test_data"
))

DamageDetective_core <- function(
    project_name,
    filtered_path,
    seurat_input   = NULL,
    min_cells      = 0,
    soupx          = FALSE,
    raw_path       = NULL,
    filter_rbc     = TRUE,
    hemo_threshold = 20,
    isolate_cd45   = FALSE,
    filter_output  = TRUE,
    output_path    = "./",
    resolution     = 0.1,
    cluster_ranks  = 1,
    organism       = "Hsap",
    verbose        = TRUE
){
  # Receive & prepare input ------------------------------------

  message("\nBeginning DamageDetective analysis for ", project_name, "...")

  # Create output directory structure
  sub_dirs <- c("RBCQC", "CellQC", "Filtered")

  for (sub_dir in sub_dirs) {
    # Generate default subdirectory names
    full_path <- file.path(output_path, sub_dir)

    # Create subdirectory if needed
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = FALSE)
    }
  }

  # Check if the input organism is valid
  if (!is.null(organism) & !organism %in% c("Hsap", "Mmus")) {
    stop("Please ensure 'organism' is one of the following 'Hsap', 'Mmus'")
  }

  # Read in required organism annotation file
  if (organism == "Hsap") {
    annotations <- get("human_annotations", envir = asNamespace("DamageDetective"))
  } else if (organism == "Mmus") {
    annotations <- get("mouse_annotations", envir = asNamespace("DamageDetective"))
  }

  # EITHER create project Seurat object with filtered counts
  if (is.null(seurat_input)) {

    # Create Seurat object using filtered counts (zipped)
    table_of_counts <- suppressWarnings(Read10X(filtered_path))

    # Create a Seurat object & record cell number
    Seurat <- suppressWarnings(CreateSeuratObject(
      counts = table_of_counts,
      project = project_name,
      min.cells = min_cells
    ))

    # Store unfiltered cell barcodes
    storage <- Seurat

    # Give the unfiltered cell number
    initial_cells <- length(Cells(Seurat))

    message("\u2714 Seurat object created")

  }

  # OR use input Seurat object if provided

  else {

   # If specified input is Seurat object, take it as such
   Seurat <- seurat_input

   # Check if the input is a Seurat object
   if (!methods::is(test_data, "Seurat")) {
     stop("Error: Ensure seurat_input is either a Seurat object or NULL.")
   }

   # Store unfiltered cell barcodes
   storage <- Seurat

   # Calculate the unfiltered cell number
   initial_cells <- length(Cells(Seurat))

   message("\u2714 Seurat object loaded")

  }



  # Optional ambient RNA correction with SoupX ------------------------------------

  if (soupx) {

    message("\u2714 Beginning SoupX correction...")

    table_of_droplets <- suppressWarnings(Read10X(raw_path))

    # Use the soupx_calculation() function to run SoupX ambient RNA correction
    Seurat <- soupx_calculation(table_of_droplets = table_of_droplets,
                                table_of_counts = table_of_counts,
                                min_cells = min_cells,
                                project_name = project_name)

    message("\u2714 SoupX correction complete")

  }

  # Default red blood cell QC ------------------------------------

  if (filter_rbc) {

    # Use the input organism for hemoglobin gene nomenclature

    if (organism == "Hsap") {
      hemo_gene <- c("HBA1", "HBA2", "HBB") # hemoglobin subunit genes
      cd45_gene <- c("PTPRC", "PTPRCAP")    # CD45 as proxy for immune cells
    }

    if (organism == "Mmus") {
      hemo_gene <- c("Hba-a1", "Hba-a2", "Hbb-bt", "Hbb-bs") # lower case & different!
      cd45_gene <- c("Ptprc", "Ptprcap")
    }

    # Calculate proportion of hemo gene expression & store as meta data
    Seurat[['hemo.percent']] <- PercentageFeatureSet(
      object   = Seurat,
      features = intersect(hemo_gene, rownames(Seurat@assays$RNA)),
      assay    = "RNA"
    )

    # Calculate proportion of CD45 gene expression & store as meta data
    Seurat[['ptprc.percent']] <- PercentageFeatureSet(
      object   = Seurat,
      features = intersect(cd45_gene, rownames(Seurat@assays$RNA)),
      assay    = "RNA"
    )

    # Label barcodes if they are likely RBCs or have hemoglobin contamination
    Seurat[['RBC']] <- ifelse(Seurat$hemo.percent >=  hemo_threshold, "RBC", "non-RBC")

    # Calculate the number of RBCs (account for none)
    RBC_number <- NULL

    if (all(Seurat$hemo.percent <= hemo_threshold)) {

      # Say there were no RBCs
      RBC_number = 0
      message("\u2714 No red blood cells found")

    } else {

      # Count RBCs if they are present
      RBC <- subset(Seurat, RBC == "RBC")
      RBC_number <- length(Cells(RBC))

      message("\u2714 Red blood cells removed")

    }

    # Create data frame with cell barcdoes and their RBC label
    rbc_df <- as.data.frame(Seurat@meta.data)
    RBC_percent <- ( RBC_number / initial_cells ) * 100

    # Plot hemoglobin vs CD45 with create_rbc_plot() function
    hemo_feature_QC <- create_rbc_plot(rbc_df = rbc_df,
                                       project_name = project_name,
                                       RBC_percent = RBC_percent,
                                       initial_cells = initial_cells)

    # Save output  to directory specified
    saveRDS(hemo_feature_QC, file.path(output_path, "RBCQC", paste0(project_name, "_RBCQC", ".rds")))

    # Filtering step after visualization, keep those NOT marked as RBC
    Seurat <- subset(Seurat, RBC == "non-RBC")

  }



  # Optional IMC QC ------------------------------------

  if (isolate_cd45) {

    # Generate IMC subdirectory name
    full_path <- file.path(output_path, "IMCQC")

    # Create IMC subdirectory if needed
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = FALSE)
    }

    # Mark cells as immune if CD45 is present in the total proporiton of genes expression (aka not zero)
    Seurat[['IMC']] <- ifelse(Seurat$ptprc.percent > 0, "IMC", "non-IMC")

    # Calculate the number of immune cells
    IMC <- subset(Seurat, IMC == "IMC")
    IMC_number <- length(Cells(IMC))

    # Plot the barcodes to visualize how many were removed (non-IMC)
    imc_df <- as.data.frame(Seurat@meta.data)

    # Calculate the percentage of total that are immune
    IMC_percent <- ( IMC_number / initial_cells ) * 100

    # Generate scatter plot with labelled IMCs using create_imc_plot() function
    immune_feature_QC <- create_imc_plot(imc_df = imc_df,
                                         project_name = project_name,
                                         IMC_percent = IMC_percent)

    # Save output to subdirectory created
    saveRDS(immune_feature_QC, file.path(output_path, "IMCQC", paste0(project_name, "_IMCQC", ".rds")))

    # Filtering step after visualization, keep those that are immune cells
    Seurat <- subset(Seurat, IMC == "IMC")

    message("\u2714 Immune cells isolated")

  }



  # Identify damaged cells with DamageDetective------------------------------------

  # Redefine starting cells to account for possible RBC & IMC filtering
  initial_cells <- length(Cells(Seurat))

  # Use the DamageDetective_calculation() function to identify damaged cells

  DamageDetective_output <- DamageDetective_calculation(organism = organism,
                                  Seurat = Seurat,
                                  resolution = resolution,
                                  cluster_ranks = cluster_ranks,
                                  annotations = annotations,
                                  initial_cells = initial_cells,
                                  project_name = project_name,
                                  output_path = output_path)

  # Add output to function environment
  Seurat  <- DamageDetective_output$Seurat
  DamageDetective <- DamageDetective_output$DamageDetective


  message("\u2714 DamageDetective damaged cell predictions")


  # Clean & save output ------------------------------------

  # Filter cells for saving the Seurat object itself (but only if specified)

  # Filter based on DamageDetective annotations
  if (filter_output) {

    # Save a list of barcodes with DamageDetective annotations (cell, damaged)
    storage_cells <- data.frame(barcode = rownames(storage@meta.data))
    clean_cells <- data.frame(barcode = rownames(Seurat@meta.data), DamageDetective = Seurat@meta.data$DamageDetective)
    write.csv(clean_cells, file = file.path(output_path, "/Filtered/", paste0(project_name, "_barcodes.csv")), row.names = FALSE, quote = FALSE)

    # Account for cells that may have been filtered
    annotated_cells <- merge(storage_cells, clean_cells, by = "barcode", all.x = TRUE)
    annotated_cells$DamageDetective[is.na(annotated_cells$DamageDetective)] <- "removed"

    write.csv(annotated_cells, file = file.path(output_path, "Filtered", paste0(project_name, "_barcodes.csv")), row.names = FALSE, quote = FALSE)


    # Filter damaged cells only according to DamageDetective estimations
    Seurat <- subset(Seurat, DamageDetective == "cell")

    # Filter unnecessary columns
    columns <- c("ptprc.percent", "RBC", "IMC", "quality", "hemo.percent",  "QC", "DamageDetective.mi", "DamageDetective.ri", "DamageDetective.complexity", "DamageDetective", "DamageDetective.droplet_qc")

    for (column in columns){
      if (column %in% colnames(Seurat@meta.data)) {
        Seurat[[column]] <- NULL }
    }

    # Save the filtered Seurat object
    saveRDS(Seurat, file.path(output_path, "Filtered", paste0(project_name, "_filtered.rds")))

  }

  # Don't filter if specified not to
  if (filter_output == FALSE) {

    # Save object without removing meta data columns (user can visualise for themsleves)

    # Save a list of barcodes with DamageDetective annotations (cell, damaged)
    storage_cells <- data.frame(barcode = rownames(storage@meta.data))
    clean_cells <- data.frame(barcode = rownames(Seurat@meta.data), DamageDetective = Seurat@meta.data$DamageDetective)
    write.csv(clean_cells, file = file.path(output_path, "/Filtered/", paste0(project_name, "_barcodes.csv")), row.names = FALSE, quote = FALSE)

    # Account for cells that may have been filtered
    annotated_cells <- merge(storage_cells, clean_cells, by = "barcode", all.x = TRUE)
    annotated_cells$DamageDetective[is.na(annotated_cells$DamageDetective)] <- "removed"

    write.csv(annotated_cells, file = file.path(output_path, "Filtered", paste0(project_name, "_barcodes.csv")), row.names = FALSE, quote = FALSE)
    saveRDS(Seurat, file.path(output_path, "Filtered", paste0(project_name, "_unfiltered.rds")))

  }

  message("\u2714 DamageDetective analysis complete\n")

  return(Seurat)

}
