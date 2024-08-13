#' limiricCore
#'
#' Helper function to enable the limiric function to have list inputs
#'
#' @name limiricCore
#'
#' @param ProjectName String with project or sample name
#' @param FilteredPath Directory of filtered alignment output
#' @param SeuratInput Seurat object to be used as input over raw files. Default NULL
#' @param MinCells In how many cells should a gene be expressed to be kept
#' @param SoupX Perform ambient RNA correction, if TRUE RawPath must be given. Default is FALSE
#' @param RawPath Directory of unfiltered alignment output
#' @param DropletQC Verify output with DropletQC, if TRUE VelocytoPath must be given. Default is FALSE
#' @param VelocytoPath Directory of Veocyto filtered alignment output
#' @param FilterRBC Whether or not red blood cells should be removed, Default is TRUE
#' @param IsolateCD45 Discard non-immune cells. Default is FALSE
#' @param FilterOutput Should output contain no damaged cells. Default is TRUE
#' @param OutputPath Directory where limiric output cen be generated
#' @param Organism "Hsap" if human sample or "Mmus" if mouse sample
#'
#' @return (list) Output list storing the final Seurat object
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
#' @export
#'
#' @examples
#' \dontrun{
#' # Detect damaged cells
#'
#' SRR1234567 <- limiricCore(
#'   ProjectName  = "SRR1234567",
#'   FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
#'   SoupX        = TRUE,
#'   RawPath      = "/home/user/alignment/SRR1234567/raw/",
#'   DropletQC    = TRUE,
#'   IsolateCD45  = TRUE,
#'   VelocytoPath = "/home/user/alignment/velocyto/",
#'   OutputPath   = "/home/user/alignment/limiric/"
#' )
#' }

utils::globalVariables(c(
  "human_annotations", "mouse_annotations", "hbb.percent", "ptprc.percent",
  "nCount_RNA", "nFeature_RNA", "gene_name", "seurat_clusters", "mt.percent",
  "rb.percent", "avg_mt_percent", "avg_rb_percent", "rank_mt", "rank_rb",
  "combined_rank", "MtRb", "quality", "QC"
))

limiricCore <- function(

  ProjectName,            # String with the name of the project
  FilteredPath,           # Path to zipped filtered output of alignment
  SeuratInput  = NULL,    # Input is a seurat object
  MinCells     = 0,       # Minimum number of cells a gene must be expressed in to be retained
  SoupX        = FALSE,   # Default is not to do SoupX ambient RNA filtering
  RawPath      = NULL,    # Path to zipped raw output of alignment
  DropletQC    = FALSE,   # Default is not to include nuclear fraction data
  VelocytoPath = NULL,    # Path to zipped filtered output of velocyto
  FilterRBC    = TRUE,    # Default is to remove red bloods cells
  IsolateCD45  = NULL,    # default is not to filter immune cells
  FilterOutput = TRUE,    # Default is to remove damaged cells
  OutputPath   = "./",    # default output directory is the working directory
  Organism     = "Hsap"   # either : human Hsap, mouse Mmus

){

  set.seed(7777)


  # Preparation ####
  # If an Organisms was added, check if is is valid
  if (!is.null(Organism) & !Organism %in% c("Hsap", "Mmus")) {
    stop("Invalid organism. Must be one of 'Hsap', 'Mmus'.")
  }



  # Read in the corresponding species gene annotation file using the specified organism (Hsap, Mmus) for correct gene nomenclature
  if (Organism == "Hsap") {
    annotations <- get("human_annotations", envir = asNamespace("limiric"))
  }

  if (Organism == "Mmus") {
    annotations <- get("mouse_annotations", envir = asNamespace("limiric"))
  }


  # Verify output directory structure & create if absent
  sub_dirs <- c("RBCQC", "CellQC", "Filtered")

  for (sub_dir in sub_dirs) {

    # Generate the names of the directories we expect
    full_path <- file.path(OutputPath, sub_dir)

    # If directory does not exist it will be created
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = FALSE)
    }

  }

  # Reading input ####

  cat("\nBeginning  limiric  analysis for", ProjectName, "...\n")

  if (!is.null(SeuratInput)) {

    # If specified input is Seurat object, take it as such
    Seurat <- SeuratInput

    storage <- Seurat # For unfiltered barcode annotations

    # Calculate the unfiltered cell number
    InitialCells <- length(Cells(Seurat))
    RhoEstimate <- NULL

    cat("\u2714 Seurat object loaded\n")
  }

  else

  {

    # Create Seurat object using filtered counts (must be zipped input files)
    TableOfCounts <- suppressWarnings(Read10X(FilteredPath))

    # Create a Seurat object & check the cell number of the sample
    Seurat <- suppressWarnings(CreateSeuratObject(
      counts = TableOfCounts,
      project = ProjectName,
      min.cells = MinCells,
      min.features = 50
    ))

    storage <- Seurat

    # Give the unfiltered cell number
    InitialCells <- length(Cells(Seurat))
    RhoEstimate <- NULL

    cat("\u2714 Seurat object created\n")
  }

  # Ambient RNA correction with SoupX ####

  if (SoupX == TRUE) {


    # Perform SoupX correction with the raw output (must be zipped)
    TableOfDroplets <- suppressWarnings(Read10X(RawPath))

    # Run Soup X
    sc <- SoupChannel(TableOfDroplets, TableOfCounts, calcSoupProfile = FALSE) # Create the soup channel (sc)

    # Estimate the contamination
    sc <- estimateSoup(sc)

    # Use Seurat to cluster the filtered matrix, although not essential it is recommended to get better estimations
    SeuratSoup <- suppressWarnings(CreateSeuratObject(TableOfCounts, min.cells = 0))
    SeuratSoup <- suppressWarnings(SCTransform(SeuratSoup, verbose = FALSE) %>%
                                     RunPCA(verbose = FALSE) %>%
                                     RunUMAP(dims = 1:30, verbose = FALSE) %>%
                                     FindNeighbors(dims = 1:30, verbose = FALSE) %>%
                                     FindClusters(verbose = FALSE))

    # Adding the cluster embeddings to the SoupX object
    MetaData <- SeuratSoup@meta.data
    UmapEmbedding <- SeuratSoup@reductions$umap@cell.embeddings
    sc <- suppressWarnings(setClusters(sc, setNames(MetaData$seurat_clusters, rownames(MetaData))))
    sc <- suppressWarnings(setDR(sc, UmapEmbedding, c("UMAP_1", "UMAP_2")))

    # With defined clusters, run the main SoupX function to calculate the contamination fraction rho where rho E (0, 1) and the closer to 1, the more contaminated
    sc <- autoEstCont(sc, verbose = FALSE, doPlot = FALSE)
    RhoEstimate <- sc$metaData$rho[1] # record the SoupX contamination estimate

    # Silencing the output ...
    # Open a connection to a temporary file for writing
    tmp_conn <- file(tempfile(), open = "wt")

    # Redirect standard output and messages to the temporary file
    sink(tmp_conn)
    sink(tmp_conn, type = "message")

    # Call the (now silenced) verbose function
    # Output integer matrix of soup-corrected reads (unzipped output) where contaminated reads are removed
    AdjMatrix <- suppressWarnings(adjustCounts(sc, roundToInt = T))

    # Reset output redirection
    sink(NULL)
    sink(NULL, type = "message")

    # Close the connection
    close(tmp_conn)

    # Closed! No more silencing (everything else is silent..)

    # Output results in Seurat object for continued workflow (will replace previous Seurat with corrected one)
    Seurat <- suppressWarnings(CreateSeuratObject(counts = AdjMatrix, # SoupX corrected count matrix
                                                  min.cells = MinCells, # At least one cell must express the gene for the gene to be included in the count matrix
                                                  min.features = 0,
                                                  project = ProjectName))

    storage <- Seurat

    cat("\u2714 SoupX correction complete\n")

  }

  # Nuclear fraction with DropletQC metric evaluation ####

  if (DropletQC == TRUE) {

    # Define parameters for file paths and gene names
    spliced_file <- file.path(VelocytoPath, "spliced.mtx.gz")
    barcodes_file <- file.path(VelocytoPath, "barcodes.tsv.gz")
    features_file <- file.path(VelocytoPath, "features.tsv.gz")
    unspliced_file <- file.path(VelocytoPath, "unspliced.mtx.gz")

    # Add data from spliced and unspliced counts
    spliced <- ReadMtx(mtx      = spliced_file,
                       cells    = barcodes_file,
                       features = features_file)

    spliced <- suppressWarnings(CreateSeuratObject(
      counts    = spliced,
      project   = "spliced",
      min.cells = MinCells))

    unspliced <- ReadMtx(mtx      = unspliced_file,
                         cells    = barcodes_file,
                         features = features_file)

    unspliced <- suppressWarnings(CreateSeuratObject(
      counts    = unspliced,
      project   = "unspliced",
      min.cells = MinCells))

    # Calculate nuclear fraction
    ExonSum   <- Matrix::colSums(spliced[['RNA']]$counts)   # summing over all the genes for each cell (1 reading per cell)
    IntronSum <- Matrix::colSums(unspliced[['RNA']]$counts)
    NuclearFraction <- IntronSum / (ExonSum + IntronSum)
    nf <- data.frame(barcode = rownames(unspliced@meta.data), nf = NuclearFraction)

    # Add nf to Seurat object (ensure row order (barcode) correct)
    Seurat@meta.data[['nf']] <- nf$nf[match(rownames(Seurat@meta.data), nf$barcode)]

    cat("\u2714 Velocyto nuclear fraction output added\n")

  }

  # Red blood cell QC ####

  if (FilterRBC == TRUE) {

    # Use the specified organism (Mmul, Hsap, Mmus) for correct gene nomenclature
    if (Organism == "Hsap") {

      # Define genes of interest
      HemoGene <- c("HBA1", "HBA2", "HBB")          # haemoglobin subunit genes
      cd45Gene <- c("PTPRC", "PTPRCAP")             # CD45 as proxy for immune cell & CD45 polypeptide-associated protein

    }


    if (Organism == "Mmus") {

      # Define genes of interest
      HemoGene <- c("Hba-a1", "Hba-a2", "Hbb-bt", "Hbb-bs") # Mmus genes are lower case & different
      cd45Gene <- c("Ptprc", "Ptprcap")

    }

    # Add metadata columns to store the expression values of the genes of interest in each barcode
    Seurat[['hbb.percent']] <- PercentageFeatureSet(
      object   = Seurat,
      features = intersect(HemoGene, rownames(Seurat@assays$RNA)),
      assay    = "RNA"
    )

    Seurat[['ptprc.percent']] <- PercentageFeatureSet(
      object   = Seurat,
      features = intersect(cd45Gene, rownames(Seurat@assays$RNA)),
      assay    = "RNA"
    )

    # RBC
    # Label barcodes if they are likely red blood cells or have red blood cell contamination
    Seurat[['RBC']] <- ifelse(Seurat$hbb.percent !=  0, "RBC", "non-RBC")

    # Calculate the number of red blood cells (account for no hemoglobin expression)
    RBCNumber <- NULL

    if (all(Seurat$hbb.percent == 0)) {

      # Say there were no RBCs
      RBCnumber = 0
      cat("\u2714 No red blood cells found\n")

    }

    else {

      # Count RBCs if they are present
      RBC <- subset(Seurat, hbb.percent > 0)
      RBCnumber <- length(Cells(RBC))

    }

    # Plot the hemoglobin percent
    rbcDf <- as.data.frame(Seurat@meta.data)
    RBC_percent <- ( RBCnumber / InitialCells ) * 100

    # Plot haemo vs CD45
    hemo_feature_QC <- ggplot(rbcDf, aes(x = hbb.percent, y = ptprc.percent, color = RBC)) +
      geom_point(size = 0.6) +
      scale_color_manual(values = c("non-RBC" = "grey","RBC" = "#5372B4")) +
      xlab("Haemoglobin expression") +
      ylab("CD45 expression") +
      labs(title = ProjectName) +
      annotate("text", x = Inf, y = Inf, label = paste("Contamination of ", round(RBC_percent, 2), "%", "for ", InitialCells, " cells"), hjust = 1.1, vjust = 2.4) +
      theme_classic() +
      theme(
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(hjust = 0.5, vjust = -1, face = "bold", size = 12),
        axis.title.y = element_text(hjust = 0.5, vjust = 3, face = "bold", size = 12),
        axis.text.x  = element_text(angle = 0, hjust = 0.5, vjust = -1.2, size = 12),
        axis.text.y  = element_text(angle = 0, hjust = 1, vjust = -1.2, size = 12),
        legend.title = element_text(face = "bold"),
        legend.position = "below",
        legend.box.background = element_rect(colour = "black"),
        strip.background = element_blank(),
        plot.margin      = unit(c(0.2, 0.2, 0.2, 0.4), 'cm'),
        plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        panel.border     = element_rect(colour = "black", fill=NA, linewidth =1)

      )

    # Save output  to directory specified (or default) : here we don't want the individual plots, we want them all together
    saveRDS(hemo_feature_QC, file.path(OutputPath, "RBCQC", paste0(ProjectName, "_RBCQC", ".rds")))

    # Filtering step after visualization, keep those not marked as RBC
    Seurat <- subset(Seurat, RBC == "non-RBC")

    cat("\u2714 Red blood cells removed\n")

  }

  # IMC QC ####

  if (IsolateCD45 == TRUE) {

    # Generate the names of the directories we expect
    full_path <- file.path(OutputPath, "IMCQC")

    # If directory does not exist it will be created
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = FALSE)
    }

    # For any numerical value, such as 0, only barcodes with CD45 percentage above this (not equal to) will be marked as IMC
    Seurat[['IMC']] <- ifelse(Seurat$ptprc.percent > 0, "IMC", "non-IMC")

    # Calculate the number of immune cells (account for no immune expression)
    IMC <- subset(Seurat, IMC == "IMC")
    IMC_number <- length(Cells(IMC))

    # Plot the barcodes to visualize how many were removed (non-IMC)
    imcDf <- as.data.frame(Seurat@meta.data)

    # Calculate the percentage of total that are immune
    IMC_percent <- ( IMC_number / InitialCells ) * 100


    # Plotting metrics that give good spread of data, UMI & gene counts per barcode
    immune_feature_QC <- ggplot(imcDf, aes(x = nCount_RNA, y = nFeature_RNA, color = IMC)) +
      geom_point(size = 0.6) +
      scale_color_manual(values = c("IMC" = "grey","non-IMC" = "#5372B4")) +
      labs(title = ProjectName) +
      annotate("text", x = Inf, y = Inf, label = paste("Percentage of immune cells: ", round(IMC_percent, 2), "%"), hjust = 1.1, vjust = 2.4) +
      theme_classic() +
      theme(
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.4), 'cm'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title.x = element_text(hjust = 0.5, vjust = -1, face = "bold", size = 12),
        axis.title.y = element_text(hjust = 0.5, vjust = 3, face = "bold", size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = -1.2, size = 12),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = -1.2, size = 12),
        legend.title = element_text(face = "bold"),
        legend.position = "below",
        legend.box.background = element_rect(colour = "black"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth =1)

      )

    # Save output  to directory specified (or default) : here we don't want the individual plots, we want them all together
    saveRDS(immune_feature_QC, file.path(OutputPath, "IMCQC", paste0(ProjectName, "_IMCQC", ".rds")))

    # Filtering step after visualization, keep those that are immune cells
    Seurat <- subset(Seurat, IMC == "IMC")

    cat("\u2714 Immune cells isolated\n")
  }

  # Identify damaged cells with limiric ####

  # Redefine starting cells to account for possible RBC & IMC filtering
  InitialCells <- length(Cells(Seurat))

  # Define mitochondrial and ribosomal low dimensional gene space

  if (Organism == "Hsap") {
    # Get gene annotations for mitochondrial (MT) & ribosomal (RB) genes
    mt_gene_annotations <- annotations[grep("MT-", annotations$gene_name, perl=TRUE),]
    mt_gene_annotations <- mt_gene_annotations[grepl("protein_coding", mt_gene_annotations$gene_biotype, perl=TRUE),]
    mt_genes <- mt_gene_annotations %>% pull(gene_name)

    # isolate ribosomal genes
    rps_genes <- annotations[grep("RPS", annotations$gene_name, perl=TRUE),]
    rps_genes <- rps_genes[grepl("protein_coding", rps_genes$gene_biotype, perl=TRUE),]
    rps_genes <- rps_genes %>% pull(gene_name)
    rpl_genes <- annotations[grep("RPL", annotations$gene_name, perl=TRUE),]
    rpl_genes <- rpl_genes[grepl("protein_coding", rpl_genes$gene_biotype, perl=TRUE),]
    rpl_genes <- rpl_genes %>% pull(gene_name)
    rb_genes  <- c(rps_genes, rpl_genes)

    # combine mt and rb genes
    mt_rb_genes <- c(mt_genes, rb_genes)
    mt_rb_genes <- unique(mt_rb_genes)
  }

  if (Organism == "Mmus") {

    # Get gene annotations for mitochondrial genes
    mt_genes <- annotations[grep("mt-", annotations$gene_name, perl = TRUE), ]
    mt_genes <- mt_genes[grepl("mitochondrially encoded", mt_genes$description, perl = TRUE), ]
    mt_genes <- mt_genes %>% pull(gene_name)


    # isolate ribosomal genes
    rb_genes <- annotations[grepl("ribosomal", annotations$description, perl = TRUE), ]
    rb_genes <- rb_genes[grepl("protein_coding", rb_genes$gene_biotype, perl = TRUE), ]
    rb_genes <- rb_genes %>% pull(gene_name)

    # combine mt and rb genes
    mt_rb_genes <- c(mt_genes, rb_genes)
    mt_rb_genes <- unique(mt_rb_genes)

  }


  # Reduce based on mt & rb genes only (this occurs in a separate Seurat object (temp))
  limiric <- subset(Seurat, features = intersect(mt_rb_genes, rownames(Seurat@assays$RNA)))

  limiric <- NormalizeData(limiric, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE) %>% # high resolution hoping to finely locate population with as little inclusion of non-damaged cells as possible
    RunTSNE(dims = 1:10, verbose = FALSE, check_duplicates = FALSE) # In this case (with far fewer genes) tSNE is prefered

  # Annotations to the reduced object based on the unreduced seurat object (all genes)
  DefaultAssay(limiric) <- "RNA"
  DefaultAssay(Seurat) <- "RNA"

  # Define mitochondrial expression in the temp Seurat object using the counts from the actual object
  limiric$mt.percent <- PercentageFeatureSet(
    object   = Seurat,
    features = intersect(mt_genes, rownames(Seurat@assays$RNA)),
    assay    = "RNA"
  )

  Seurat$limiric.mi <- limiric$mt.percent # Transfer to Seurat

  # Define ribosomal expression
  limiric$rb.percent <- PercentageFeatureSet(
    object   = Seurat,
    features = intersect(rb_genes, rownames(Seurat@assays$RNA)),
    assay    = "RNA"
  )

  Seurat$limiric.ri <- limiric$rb.percent # Transfer to Seurat

  # Calculate complexity
  matrix <- GetAssayData(limiric, layer = "data")
  complexity_metric <- colSums(matrix > 0)
  limiric$complexity <- complexity_metric
  Seurat$limiric.complexity <- limiric$complexity # Transfer to Seurat



  # Automatically find the damaged cell population
  # Calculate the average mt.percent and rb.percent for each cluster, add ranks, and identify the best cluster
  best_cluster <- limiric@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      avg_mt_percent = mean(mt.percent, na.rm = TRUE),
      avg_rb_percent = mean(rb.percent, na.rm = TRUE)
    ) %>%
    mutate(
      rank_mt = rank(-avg_mt_percent),  # Negative for descending order (highest mt to be #1)
      rank_rb = rank(avg_rb_percent),
      combined_rank = rank_mt + rank_rb
    ) %>%
    arrange(combined_rank) %>%
    slice(1) %>%
    pull(seurat_clusters) %>%
    as.character()

  # Label all cells belonging to this cluster as "damaged
  limiric$seurat_clusters <- ifelse(limiric$seurat_clusters == best_cluster,'damaged', limiric$seurat_clusters)

  # Add Cell QC meta data to object
  # MtRb results
  limiric$MtRb <- ifelse(limiric$seurat_clusters == "damaged", "damaged", "cell")

  # Calculated percentage of damaged cells
  damaged <- subset(limiric, MtRb == "damaged")
  damagedCells <- length(Cells(damaged))
  damagedPercent <- (damagedCells / InitialCells) * 100

  # Add cluster #s to actual seurat object
  Seurat$limiric <- limiric$MtRb



  # Visualise limiric output
  # tSNE of mt percent
  mt_plot <- FeaturePlot(limiric, features = c("mt.percent"), cols = c("#E1E1E1", "#5372B4"), pt.size = 1) +
    NoAxes() + labs(caption = "Mitochondrial gene expression") +
    theme(
      plot.title   = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))

  # tSNE of rb percent
  rb_plot <- FeaturePlot(limiric, features = c("rb.percent"), cols = c("#E1E1E1", "#5372B4"), pt.size = 1) +
    NoAxes() + labs(caption = "Ribosomal gene expression") +
    theme(
      plot.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))


  # View labels in reduced space
  cluster_plot <- DimPlot(
    limiric, pt.size = 1, group.by = "MtRb", cols = c("cell" = "grey", "damaged" = "#5372B4")) +
    labs(caption = expression("Damaged cells identified by " * italic("limiric"))) + NoAxes() +
    theme(
      plot.title = element_blank(),
      plot.caption = element_text(hjust = 0.5, size = 16),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1)
    )

  # Complexity plot
  # tSNE of mt percent
  complexity_plot <- FeaturePlot(limiric, features = c("complexity"), cols = c("#E1E1E1", "#5372B4"), pt.size = 1) +
    NoAxes() + labs(caption = "Complexity score") +
    theme(
      plot.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))


  # FANCY PLOT COMBINING HERE

  title <- ProjectName
  label <- paste("Estimated damaged cells: ", round(damagedPercent, 2), "%, ", InitialCells, " cells")
  dim <- readPNG(system.file("extdata", "tSNE.png", package = "limiric"))

  # Combine the plots
  limiric_plot <- plot_grid(mt_plot, complexity_plot, rb_plot, cluster_plot, ncol = 2)
  limiric_plot <- (mt_plot | complexity_plot) / (rb_plot | cluster_plot)

  # Create the title and subtitle
  title <- ggdraw() + draw_label(ProjectName, fontface = 'bold', x = 0.45, y = 0.1, hjust = 0.5, size = 20)
  subtitle <- ggdraw() + draw_label(paste("Estimated", round(damagedPercent, 2), "% damaged of ", InitialCells, " cells"), x = 0.45, hjust = 0.5, size = 16)

  # Create the logo and dim plots
  logo_plot <- ggdraw() + draw_image(logo, x = 0.25, y = 0.3, width = 0.4, height = 0.4)
  dim_plot <- ggdraw() + draw_image(dim, x = 0.54, y = 0.5, width = 0.8, height = 0.8)

  # Combine everything into the final plot
  final_plot <- plot_grid(
    title, subtitle, limiric_plot, dim_plot,
    ncol = 1,
    rel_heights = c(0.1, 0.1, 1, 0.2)
  )

  # Set the background to white for the entire plot
  final_plot <- final_plot + theme(plot.background = element_rect(fill = "white", color = "white"))

  # Save the final plot
  ggsave(file.path(OutputPath, "/CellQC/", paste0(ProjectName, ".png")), plot = final_plot, width = 12, height = 10, dpi = 300)

  cat("\u2714 limiric  damaged cell predictions\n")

  # Combine with DropletQC findings ####

  # Only continue if velocyto present
  if (DropletQC == TRUE) {

    # Generate the names of the directories we expect
    full_path <- file.path(OutputPath, "DropletQC")

    # If directory does not exist it will be created
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = FALSE)
    }


    # Identify droplets to be removed with DropletQC
    # Extract nf meta data & associated cell barcode from Seurat object
    edDf <- data.frame(nf = as.numeric(Seurat$nf), umi = Seurat$nCount_RNA)

    # Use DropletQC function to identify empty droplets
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

    # Add DropletQC output to Seurat object metadata
    Seurat$DropletQC <- dcresultsDf$df$cell_status

    # Combine limiric and DropletQC output
    # Initialise an empty metadata column
    Seurat$QC <- "cell" # all barcodes marked as cell by default at the start

    Seurat$QC <- Seurat@meta.data %>%
      dplyr::mutate(QC = case_when(
        DropletQC == "cell" & limiric == "cell" ~ "cell",
        DropletQC == "damaged_cell" & limiric == "damaged" ~ "agreed",
        DropletQC == "damaged_cell" & limiric == "cell" ~ "DropletQC",
        DropletQC == "cell" & limiric == "damaged" ~ "limiric",
        DropletQC == "empty_droplet" ~ "cell" # For our case, only interested in damaged cells not empty droplets
      ))

    # Transfer labels to reduced object to visualise labels
    limiric$QC <- Seurat$QC

    # Label damaged cells & empty drops
    # Unfiltered_cells <- length(Cells(Seurat))

    Seurat$quality <- ifelse(Seurat$QC == "agreed", "discard", "retain") # only discarding those marked as damaged by both metrics

    # Quantify damaged cells
    length_total <- length(Cells(Seurat))
    length_disagree <- subset(Seurat, quality == "retain")
    length_disagree <- length(Cells(length_disagree))

    if (length_total == length_disagree) {

      length_DropletQC <- subset(Seurat, DropletQC == "cell" | DropletQC == "empty_droplet")
      length_DropletQC <- length(Cells(length_DropletQC))
      length_damaged_DropletQC <- (length_total - length_DropletQC)

      length_limiric <- subset(Seurat, limiric == "cell")
      length_limiric <- length(Cells(length_limiric))
      length_damaged_limiric <- (length_total - length_limiric)


      DamagedPercent <- 0

      cat(paste0("\u2714 No agreement detected between ", length_damaged_DropletQC, "  DropletQC and ", length_damaged_limiric, "  limiric  damaged cells\n"))

    }

    else

    {

      seurat_damaged <- subset(Seurat, quality == "discard")
      damaged <- length(Cells(seurat_damaged))
      DamagedPercent <- ( damaged / InitialCells ) * 100

      cat("\u2714  DropletQC agreement calculated\n")

    }

    # VISUALISE

    # View labels in reduced space
    cluster_plot <- DimPlot(
      limiric, pt.size = 1, group.by = "QC", cols = c("cell" = "#D5D5D5", "agreed" = "#506BB0", "DropletQC" = "#D98C25", "limiric" = "#A4BAF2")) +
      NoAxes() +
      theme(
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth =1)
      )

    # Save the final plot
    ggsave(file.path(OutputPath, "DropletQC", paste0(ProjectName, "_clusters.png")), plot = cluster_plot, width = 5, height = 3, dpi = 300)


    # View labels in scatter plot nf vs UMI_count (DropletQC metric)
    df <- Seurat@meta.data
    df$nf <- as.numeric(df$nf)
    df$nCount_RNA <- as.numeric(df$nCount_RNA)

    scatter <- ggplot(df, aes(x = nf, y = nCount_RNA, color = QC)) +
      geom_point(size = 1.2) +
      scale_color_manual(values = c("cell" = "#D5D5D5", "agreed" = "#506BB0", "DropletQC" = "#D98C25", "limiric" = "#A4BAF2")) +
      coord_trans(y = "log10") +
      xlab("nf") + ylab("log10(UMI)") + #labs(caption = "Agreement between **_DropletQC_** and *limiric* estimations") +
      labs(caption = expression("Agreement between " * italic("DropletQC") * " and " * italic("limiric") * " estimations")) +
      theme(
        plot.caption = element_text(hjust = 0.5, vjust = -0.5, size = 16),
        axis.text = element_text(size = 16),
        axis.title.x = element_text(hjust = 0.5, vjust = -0.2, size = 16),
        axis.title.y = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "below",
        legend.box.background = element_rect(colour = "black"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
        axis.line = element_line(colour = "black")
      )

    # FANCY PLOT COMBINING
    clusters <- readPNG(file.path(OutputPath, "DropletQC", paste0(ProjectName, "_clusters.png")))

    # Create the title and subtitle
    title <- ggdraw() + draw_label(ProjectName, fontface = 'bold', x = 0.5, y = 0.1, hjust = 0.5, size = 20)
    subtitle <- ggdraw() + draw_label(paste("Agreed estimated", round(DamagedPercent, 2), "% damaged of ", InitialCells, "cells"), x = 0.5, hjust = 0.5, size = 16)

    # Create the dim plots
    cluster_plot <- ggdraw() + draw_image(clusters, x = -0.4, y = 5, width = 2.4, height = 2.4)

    # Combine everything into the final plot
    final_plot <- plot_grid(
      title, subtitle, scatter, cluster_plot,
      ncol = 1,
      rel_heights = c(0.1, 0.1, 1, 0.15)
    )

    # Set the background to white for the entire plot
    final_plot <- final_plot + theme(plot.background = element_rect(fill = "white", color = "white"))

    # Save the final plot
    ggsave(file.path(OutputPath, "DropletQC", paste0(ProjectName, ".png")), plot = final_plot, width = 12, height = 10, dpi = 300)


    # Finalise output
    # Save a list of barcodes with QC annnotations (cell, agreed, DropletQC, limiric)

    storageCells <- data.frame(barcode = rownames(storage@meta.data))
    cleanCells <- data.frame(barcode = rownames(limiric@meta.data), QCannotation = limiric@meta.data$QC)

    # Account for cells that may have been filtered
    annotatedCells <- merge(storageCells, cleanCells, by = "barcode", all.x = TRUE)
    annotatedCells$QCannotation[is.na(annotatedCells$QCannotation)] <- "removed"

    write.csv(annotatedCells, file = file.path(OutputPath, "DropletQC", paste0(ProjectName, "_barcodes.csv")), row.names = FALSE)

    # Rename column
    Seurat$limiric.DropletQC <- Seurat$QC

  }

  #  Clean & save output #####


  # Filter cells for saving the Seurat object itself (but only if specified)

  # Filter based on agreement
  if (FilterOutput == TRUE & DropletQC == TRUE) {

    # Save a list of barcodes with limiric annotations (cell, damaged)
    cleanCells <- data.frame(barcode = rownames(Seurat@meta.data), limiric = Seurat@meta.data$limiric)
    write.csv(cleanCells, file = file.path(OutputPath, "Filtered", paste0(ProjectName, "_barcodes.csv")), row.names = FALSE)

    Seurat <- subset(Seurat, QC != "agreed")

    # Filter unnecessary columns
    columns <- c("nf", "ptprc.percent", "RBC", "IMC", "quality", "nf", "hbb.percent",  "QC", "DropletQC", "limiric.mi", "limiric.ri", "limiric.complexity", "limiric", "limiric.DropletQC")

    for (column in columns){
      if (column %in% colnames(Seurat@meta.data)) {
        Seurat[[column]] <- NULL }
    }


    # Save the filtered Seurat object
    saveRDS(Seurat, file.path(OutputPath, "Filtered", paste0(ProjectName, "_agreement_filtered.rds")))

  }

  # Filter based only on limiric
  if (FilterOutput == TRUE & DropletQC != TRUE) {

    # Save a list of barcodes with limiric annotations (cell, damaged)
    storageCells <- data.frame(barcode = rownames(storage@meta.data))
    cleanCells <- data.frame(barcode = rownames(Seurat@meta.data), limiric = Seurat@meta.data$limiric)
    write.csv(cleanCells, file = file.path(OutputPath, "/Filtered/", paste0(ProjectName, "_barcodes.csv")), row.names = FALSE)

    # Account for cells that may have been filtered
    annotatedCells <- merge(storageCells, cleanCells, by = "barcode", all.x = TRUE)
    annotatedCells$limiric[is.na(annotatedCells$limiric)] <- "removed"

    write.csv(annotatedCells, file = file.path(OutputPath, "Filtered", paste0(ProjectName, "_barcodes.csv")), row.names = FALSE)


    # Filter damaged cells only according to limiric estimations
    Seurat <- subset(Seurat, limiric == "cell")

    # Filter unnecessary columns
    columns <- c("nf", "ptprc.percent", "RBC", "IMC", "quality", "nf", "hbb.percent",  "QC", "DropletQC", "limiric.mi", "limiric.ri", "limiric.complexity", "limiric", "limiric.DropletQC")

    for (column in columns){
      if (column %in% colnames(Seurat@meta.data)) {
        Seurat[[column]] <- NULL }
    }


    # Save the filtered Seurat object
    saveRDS(Seurat, file.path(OutputPath, "Filtered", paste0(ProjectName, "_filtered.rds")))

  }

  if (FilterOutput == FALSE) {

    # Save object without removing meta data columns (user can visualise for themsleves)

    # Save a list of barcodes with limiric annotations (cell, damaged)
    storageCells <- data.frame(barcode = rownames(storage@meta.data))
    cleanCells <- data.frame(barcode = rownames(Seurat@meta.data), limiric = Seurat@meta.data$limiric)
    write.csv(cleanCells, file = file.path(OutputPath, "/Filtered/", paste0(ProjectName, "_barcodes.csv")), row.names = FALSE)

    # Account for cells that may have been filtered
    annotatedCells <- merge(storageCells, cleanCells, by = "barcode", all.x = TRUE)
    annotatedCells$limiric[is.na(annotatedCells$limiric)] <- "removed"

    write.csv(annotatedCells, file = file.path(OutputPath, "Filtered", paste0(ProjectName, "_barcodes.csv")), row.names = FALSE)
    saveRDS(Seurat, file.path(OutputPath, "Filtered", paste0(ProjectName, "_unfiltered.rds")))

  }

  cat("\u2714 limiric  analysis complete.\n\n")

  return(Seurat)

}
