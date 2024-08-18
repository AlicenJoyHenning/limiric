#' limiric_core
#'
#' @name limiric_core
#'
#' @description Core helper function enabling the limiric modules to run
#'
#' @param project_name String with project or sample name
#' @param filtered_path Directory of filtered alignment output
#' @param seurat_input Seurat object to be used as input over raw files. Default NULL
#' @param min_cells In how many cells should a gene be expressed to be kept
#' @param soupx Perform ambient RNA correction, if TRUE raw_path must be given. Default is FALSE
#' @param raw_path Directory of unfiltered alignment output
#' @param droplet_qc Verify output with droplet_qc, if TRUE velocyto_path must be given. Default is FALSE
#' @param velocyto_path Directory of Veocyto filtered alignment output
#' @param filter_rbc Whether or not red blood cells should be removed, Default is TRUE
#' @param isolate_cd45 Discard non-immune cells. Default is FALSE
#' @param filter_output Should output contain no damaged cells. Default is TRUE
#' @param output_path Directory where limiric output cen be generated
#' @param organism "Hsap" if human sample or "Mmus" if mouse sample
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
#' \dontrun{
#'
#' # Detect damaged cells
#'
#' SRR1234567 <- limiric_core(
#'   project_name  = "SRR1234567",
#'   filtered_path = "/home/user/alignment/SRR1234567/filtered/",
#'   soupx        = TRUE,
#'   raw_path      = "/home/user/alignment/SRR1234567/raw/",
#'   droplet_qc    = TRUE,
#'   isolate_cd45  = TRUE,
#'   velocyto_path = "/home/user/alignment/velocyto/",
#'   output_path   = "/home/user/alignment/limiric/"
#' )
#'
#' }
#'
#' @keywords internal

utils::globalVariables(c(
  "human_annotations", "mouse_annotations", "hemo.percent", "ptprc.percent",
  "nCount_RNA", "nFeature_RNA", "gene_name", "seurat_clusters", "mt.percent",
  "rb.percent", "avg_mt_percent", "avg_rb_percent", "rank_mt", "rank_rb",
  "combined_rank", "MtRb", "quality", "QC"
))

limiric_core <- function(
    project_name,
    filtered_path,
    seurat_input  = NULL,
    min_cells     = 0,
    soupx         = FALSE,
    raw_path      = NULL,
    droplet_qc    = FALSE,
    velocyto_path = NULL,
    filter_rbc    = TRUE,
    isolate_cd45  = FALSE,
    filter_output = TRUE,
    output_path   = "./",
    organism      = "Hsap"
){
  # Receive & prepare input ------------------------------------

  cat("\nBeginning  limiric  analysis for", project_name, "...\n")

  # Ensure calculations are reproducible
  set.seed(7777)

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
    annotations <- get("human_annotations", envir = asNamespace("limiric"))
  } else if (organism == "Mmus") {
    annotations <- get("mouse_annotations", envir = asNamespace("limiric"))
  }

  # EITHER create project Seurat object with filtered counts
  if (is.null(seurat_input)) {

    # Create Seurat object using filtered counts (zipped)
    table_of_counts <- suppressWarnings(Read10X(filtered_path))

    # Create a Seurat object & record cell number
    Seurat <- suppressWarnings(CreateSeuratObject(
      counts = table_of_counts,
      project = project_name,
      min.cells = min_cells,
      min.features = 50
    ))

    # Store unfiltered cell barcodes
    storage <- Seurat

    # Give the unfiltered cell number
    initial_cells <- length(Cells(Seurat))

    cat("\u2714 Seurat object created\n")

  }

  # OR use input Seurat object if provided

  else {
   # If specified input is Seurat object, take it as such
   Seurat <- seurat_input

   # Store unfiltered cell barcodes
   storage <- Seurat

   # Calculate the unfiltered cell number
   initial_cells <- length(Cells(Seurat))

   cat("\u2714 Seurat object loaded\n")

  }



  # Optional ambient RNA correction with SoupX ------------------------------------

  if (soupx) {

     cat("\u2714 Beginning SoupX correction...\n")
    
    # Use the soupx_calculation() function to run SoupX ambient RNA correction
    Seurat <- soupx_calculation(raw_path = raw_path,
                                table_of_counts = table_of_counts,
                                min_cells = min_cells,
                                project_name = project_name)

    cat("\u2714 SoupX correction complete\n")

  }



  # Optional add nuclear fraction data for downstream DropletQC ------------------------------------

  if (droplet_qc) {

    # Define parameters for file paths and gene names
    spliced_file   <- file.path(velocyto_path, "spliced.mtx.gz")
    barcodes_file  <- file.path(velocyto_path, "barcodes.tsv.gz")
    features_file  <- file.path(velocyto_path, "features.tsv.gz")
    unspliced_file <- file.path(velocyto_path, "unspliced.mtx.gz")

    # Add data from spliced and unspliced counts to make Seurat objects
    spliced <- ReadMtx(mtx      = spliced_file,
                       cells    = barcodes_file,
                       features = features_file)

    spliced <- suppressWarnings(CreateSeuratObject(
      counts    = spliced,
      project   = "spliced",
      min.cells = min_cells))

    unspliced <- ReadMtx(mtx      = unspliced_file,
                         cells    = barcodes_file,
                         features = features_file)

    unspliced <- suppressWarnings(CreateSeuratObject(
      counts    = unspliced,
      project   = "unspliced",
      min.cells = min_cells))

    # Calculate nuclear fraction
    ExonSum   <- Matrix::colSums(spliced[['RNA']]$counts)   # summing over all the genes for each cell (1 reading per cell)
    IntronSum <- Matrix::colSums(unspliced[['RNA']]$counts)
    NuclearFraction <- IntronSum / (ExonSum + IntronSum)
    nf <- data.frame(barcode = rownames(unspliced@meta.data), nf = NuclearFraction)

    # Add nf to Seurat object (ensure row order (barcode) correct)
    Seurat@meta.data[['nf']] <- nf$nf[match(rownames(Seurat@meta.data), nf$barcode)]

    cat("\u2714 Velocyto nuclear fraction output added\n")

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
    Seurat[['RBC']] <- ifelse(Seurat$hemo.percent !=  0, "RBC", "non-RBC")

    # Calculate the number of RBCs (account for none)
    RBC_number <- NULL

    if (all(Seurat$hemo.percent == 0)) {

      # Say there were no RBCs
      RBC_number = 0
      cat("\u2714 No red blood cells found\n")

    } else {

      # Count RBCs if they are present
      RBC <- subset(Seurat, RBC == "RBC")
      RBC_number <- length(Cells(RBC))

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

    cat("\u2714 Red blood cells removed\n")

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

    cat("\u2714 Immune cells isolated\n")

  }



  # Identify damaged cells with limiric ------------------------------------

  # Redefine starting cells to account for possible RBC & IMC filtering
  initial_cells <- length(Cells(Seurat))

  # Use the limiric_calculation() function to identify damaged cells
  cat("\n Identifying damaged cells")
  limiric_output <- limiric_calculation(organism = organism,
                                  Seurat = Seurat,
                                  annotations = annotations,
                                  initial_cells = initial_cells,
                                  project_name = project_name,
                                  output_path = output_path)

  # Add output to function environment
  Seurat <- limiric_output$Seurat
  limiric <- limiric_output$limiric


  cat("\u2714 limiric  damaged cell predictions\n")



  # Combine with damaged cells identified with DropletQC ------------------------------------

  # Perform DropletQC calculations if specified
  if (droplet_qc == TRUE) {

    # Create Droplet QC subdirectories if needed
    full_path <- file.path(output_path, "droplet_qc")

    # If subdirectory does not exist it will be created
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = FALSE)
    }

    # Identify droplets to be removed with DropletQC
    Seurat <- dropletqc_calculation(Seurat)

    # Combine limiric and DropletQC output
    Seurat$QC <- "cell" # all barcodes marked as cell by default at the start

    Seurat$QC <- Seurat@meta.data %>%
      dplyr::mutate(QC = case_when(
        droplet_qc == "cell" & limiric == "cell" ~ "cell",
        droplet_qc == "damaged_cell" & limiric == "damaged" ~ "agreed",
        droplet_qc == "damaged_cell" & limiric == "cell" ~ "droplet_qc",
        droplet_qc == "cell" & limiric == "damaged" ~ "limiric",
        droplet_qc == "empty_droplet" ~ "cell" # For our case, only interested in damaged cells not empty droplets
      ))

    # Transfer labels to reduced object to visualise
    limiric$QC <- Seurat$QC
    Seurat$quality <- ifelse(Seurat$QC == "agreed", "discard", "retain") # only discarding those marked as damaged by both metrics

    # Quantify damaged cells
    length_total <- length(Cells(Seurat))
    length_disagree <- subset(Seurat, quality == "retain")
    length_disagree <- length(Cells(length_disagree))

    # Accounting for no overlap
    if (length_total == length_disagree) {

      length_droplet_qc <- subset(Seurat, droplet_qc == "cell" | droplet_qc == "empty_droplet")
      length_droplet_qc <- length(Cells(length_droplet_qc))
      length_damaged_droplet_qc <- (length_total - length_droplet_qc)

      length_limiric <- subset(Seurat, limiric == "cell")
      length_limiric <- length(Cells(length_limiric))
      length_damaged_limiric <- (length_total - length_limiric)

      damaged_percent <- 0

      cat(paste0("\u2714 No agreement detected between ", length_damaged_droplet_qc, "  droplet_qc and ", length_damaged_limiric, "  limiric  damaged cells\n"))

    }

    else {

      seurat_damaged <- subset(Seurat, quality == "discard")
      damaged <- length(Cells(seurat_damaged))
      damaged_percent <- ( damaged / initial_cells ) * 100

      cat("\u2714  DropletQC agreement calculated\n")

    }

    # Visualise agreement output
    create_dropletqc_plot(Seurat = Seurat,
                          limiric = limiric,
                          output_path = output_path,
                          project_name = project_name,
                          damaged_percent = damaged_percent,
                          initial_cells = initial_cells)


    # Save a list of barcodes with QC annnotations (cell, agreed, droplet_qc, limiric)
    storage_cells <- data.frame(barcode = rownames(storage@meta.data))
    clean_cells <- data.frame(barcode = rownames(limiric@meta.data), QC_annotation = limiric@meta.data$QC)

    # Account for cells that may have been filtered
    annotated_cells <- merge(storage_cells, clean_cells, by = "barcode", all.x = TRUE)
    annotated_cells$QC_annotation[is.na(annotated_cells$QC_annotation)] <- "removed"

    write.csv(annotated_cells, file = file.path(output_path, "droplet_qc", paste0(project_name, "_barcodes.csv")), row.names = FALSE)

    # Rename column
    Seurat$limiric.droplet_qc <- Seurat$QC

  }



  # Clean & save output ------------------------------------

  # Filter cells for saving the Seurat object itself (but only if specified)

  # If DropletQC agreement specified, filter accordingly
  if (filter_output == TRUE & droplet_qc == TRUE) {

    # Save a list of barcodes with limiric annotations (cell, damaged)
    clean_cells <- data.frame(barcode = rownames(Seurat@meta.data), limiric = Seurat@meta.data$limiric)
    write.csv(clean_cells, file = file.path(output_path, "Filtered", paste0(project_name, "_barcodes.csv")), row.names = FALSE)

    Seurat <- subset(Seurat, QC != "agreed")

    # Filter unnecessary columns
    columns <- c("nf", "ptprc.percent", "RBC", "IMC", "quality", "nf", "hemo.percent",  "QC", "droplet_qc", "limiric.mi", "limiric.ri", "limiric.complexity", "limiric", "limiric.droplet_qc")

    for (column in columns){
      if (column %in% colnames(Seurat@meta.data)) {
        Seurat[[column]] <- NULL }
    }

    # Save the filtered Seurat object
    saveRDS(Seurat, file.path(output_path, "Filtered", paste0(project_name, "_agreement_filtered.rds")))

  }

  # If just limiric, filter only on limiric annotations
  if (filter_output == TRUE & droplet_qc != TRUE) {

    # Save a list of barcodes with limiric annotations (cell, damaged)
    storage_cells <- data.frame(barcode = rownames(storage@meta.data))
    clean_cells <- data.frame(barcode = rownames(Seurat@meta.data), limiric = Seurat@meta.data$limiric)
    write.csv(clean_cells, file = file.path(output_path, "/Filtered/", paste0(project_name, "_barcodes.csv")), row.names = FALSE)

    # Account for cells that may have been filtered
    annotated_cells <- merge(storage_cells, clean_cells, by = "barcode", all.x = TRUE)
    annotated_cells$limiric[is.na(annotated_cells$limiric)] <- "removed"

    write.csv(annotated_cells, file = file.path(output_path, "Filtered", paste0(project_name, "_barcodes.csv")), row.names = FALSE)


    # Filter damaged cells only according to limiric estimations
    Seurat <- subset(Seurat, limiric == "cell")

    # Filter unnecessary columns
    columns <- c("nf", "ptprc.percent", "RBC", "IMC", "quality", "nf", "hemo.percent",  "QC", "droplet_qc", "limiric.mi", "limiric.ri", "limiric.complexity", "limiric", "limiric.droplet_qc")

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

    # Save a list of barcodes with limiric annotations (cell, damaged)
    storage_cells <- data.frame(barcode = rownames(storage@meta.data))
    clean_cells <- data.frame(barcode = rownames(Seurat@meta.data), limiric = Seurat@meta.data$limiric)
    write.csv(clean_cells, file = file.path(output_path, "/Filtered/", paste0(project_name, "_barcodes.csv")), row.names = FALSE)

    # Account for cells that may have been filtered
    annotated_cells <- merge(storage_cells, clean_cells, by = "barcode", all.x = TRUE)
    annotated_cells$limiric[is.na(annotated_cells$limiric)] <- "removed"

    write.csv(annotated_cells, file = file.path(output_path, "Filtered", paste0(project_name, "_barcodes.csv")), row.names = FALSE)
    saveRDS(Seurat, file.path(output_path, "Filtered", paste0(project_name, "_unfiltered.rds")))

  }

  cat("\u2714 limiric  analysis complete.\n\n")

  return(Seurat)

}
