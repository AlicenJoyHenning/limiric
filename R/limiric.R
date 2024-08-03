#' limiric
#'
#' @description Function for scRNA-seq damaged cell detection through low dimension mitochondrial and ribosomal cluster selection
#'
#' @param ProjectName (character) Name of the project or sample being analysed, for example "SRR1234567"
#' @param FilteredPath (character) Path to the directory where gzipped, filtered alignment output is stored, should contain matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz. For example, "/home/user/alignmentoutput/SRR1234567/filtered/"
#' @param SeuratInput (logical) Whether or not the input is a Seurat object
#' @param MinCells (numeric) When reading in the scRNA-seq data, how many cells should a gene be expressed in to be retained in the count matrix? For example, 1
#' @param SoupX (logical) Whether or not ambient RNA correction should be performed, requiring RawPath input
#' @param RawPath (character) Path to the directory where gzipped, raw alignment output is stored, should contain matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz. For example, "/home/user/alignmentoutput/SRR1234567/raw/"
#' @param DropletQC (logical) Whether or not damaged cell annotation should be compared to DropletQC damaged cell annotations, requiring VelocytoPath
#' @param VelocytoPath (character) Path to the directory where gzipped, filtered output from Velocyto is stored, should contain spliced.mtx.gz, unspliced.mtx.gz, features.tsv.gz, and barcodes.tsv.gz. For example, "/home/user/alignmentoutput/SRR1234567/velocyto/"
#' @param FilterRBC (logical) Whether or not red blood cells should be removed
#' @param IsolateCD45 (logical) Whether or not immune cells should be isolated
#' @param FilterOutput (logical) Whether or not output should filter out damaged cells, or keep all cells with their limiric annotation
#' @param OutputPath (character) Path to the directory where limiric output can be stored
#' @param Organism (character) Either "Hsap" if human sample or "Mmus" if mouse
#' @param sample_list (logical) Whether or not the input is for multiple samples defined in a list
#'
#' @return (list) Output list storing the final Seurat object
#'
#' @examples
#'
#' # Detect damaged cells in one sample after ambient RNA correction, filtering red blood cells, and isolating immune cells
#' # The sample is of human origin (default = "Hsap") and red blood cells need to be removed (default = TRUE)
#' SRR1234567 <- limiric(
#'  ProjectName  = "SRR1234567",
#'  FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
#'  SoupX        = TRUE,
#'  RawPath      = "/home/user/alignment/SRR1234567/raw/",
#'  DropletQC    = TRUE,
#'  IsolateCD45  = TRUE,
#'  VelocytoPath = "/home/user/alignment/velocyto/",
#'  OutputPath   = "/home/user/alignment/limiric/"
#' )
#'
#' # Conditions as above but for multiple samples
#' # Define samples in a list
#'
#' sample_list <- list(
#'
#' list(ProjectName  = "SRR1234567",
#'      FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
#'      SoupX        = TRUE,
#'      RawPath      = "/home/user/alignment/SRR1234567/raw/",
#'      DropletQC    = TRUE,
#'      IsolateCD45  = TRUE,
#'      VelocytoPath = "/home/user/alignment/velocyto/",
#'      OutputPath   = "/home/user/alignment/limiric/"
#'     ),
#'
#' list(ProjectName  = "SRR1234568",
#'      FilteredPath = "/home/user/alignment/SRR1234568/filtered/",
#'      SoupX        = TRUE,
#'      RawPath      = "/home/user/alignment/SRR1234568/raw/",
#'      DropletQC    = TRUE,
#'      IsolateCD45  = TRUE,
#'      VelocytoPath = "/home/user/alignment/velocyto/",
#'      OutputPath   = "/home/user/alignment/limiric/"
#'     ),
#'
#' list(ProjectName  = "SRR1234569",
#'      FilteredPath = "/home/user/alignment/SRR1234569/filtered/",
#'      SoupX        = TRUE,
#'      RawPath      = "/home/user/alignment/SRR1234569/raw/",
#'      DropletQC    = TRUE,
#'      IsolateCD45  = TRUE,
#'      VelocytoPath = "/home/user/alignment/velocyto/",
#'      OutputPath   = "/home/user/alignment/limiric/"
#'     )
#'  )
#'
#'  @import Seurat
#'  @import Matrix
#'  @import data.table
#'  @import SoupX
#'  @import DropletQC
#'  @import dplyr
#'  @import ggplot2
#'  @import cowplot
#'  @import png
#'
#' @export


# Install packages ####

# List of required packages
# required_packages <- c(
#   "Matrix", "data.table", "Seurat", "SoupX", "DropletQC",
#   "dplyr", "ggplot2", "cowplot", "png"
# )
#
# # Function to check and install missing packages
# for (pkg in required_packages) {
#     if (!require(pkg, character.only = TRUE)) {
#       install.packages(pkg, dependencies = TRUE)
#       library(pkg, character.only = TRUE)
#     }
# }
#
#
# library(Matrix)        # Helper nuclear fraction function
# library(data.table)
# library(Seurat)
# library(SoupX)         # Ambient RNA correction
# library(Seurat)        # All round backbone of analysis
# library(DropletQC)     # Identifying low quality cells (empty & damaged)
# library(dplyr)         # To use multiple functions together %>%
# library(ggplot2)       # To create "diagnostic" plots
# library(cowplot)       # To make plot grids
# library(png)



# Core function ####

limiricCore <- function(

  ProjectName,            # String with the name of the project
  FilteredPath,           # Path to zipped filtered output of alignment
  SeuratInput  = FALSE,   # Whether or not the input is a seurat object
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

  # Preparation ####
  # If an Organisms was added, check if is is valid
  if (!is.null(Organism) & !Organism %in% c("Hsap", "Mmus")) {
    stop("Invalid organism. Must be one of 'Hsap', 'Mmus'.")
  }

  # Read in the corresponding species gene annotation file using the specified organism (Hsap, Mmus) for correct gene nomenclature
  if (Organism == "Hsap") {annotations <- readRDS("/home/alicen/R/human_annotations.rds")}
  if (Organism == "Mmus") {annotations <- readRDS("/home/alicen/R/mouse_annotations.rds")}


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

  cat("\nBeginning \033[3mlimiric\033[23m analysis for", ProjectName, "...\n")

  if (SeuratInput == TRUE) {

    # If specified input is Seurat object, take it as such
    Seurat <- SeuratInput

    # Calculate the unfiltered cell number
    InitialCells <- length(Cells(Seurat))
    RhoEstimate <- NULL

    cat("✔ \033[1mSeurat\033[22m object loaded\n")
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

    # Give the unfiltered cell number
    InitialCells <- length(Cells(Seurat))
    RhoEstimate <- NULL

    cat("✔  \033[1mSeurat\033[22m object created\n")
  }

  # Ambient RNA correction with SoupX ####

  if (SoupX == TRUE) {

    cat("✔  Performing \033[1mSoupX\033[22m ambient RNA correction\n")

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

    cat("✔  \033[1mSoupX\033[22m correction complete\n")

  }

  # Nuclear fraction with DropletQC metric evaluation ####

  if (DropletQC == TRUE) {

    # Add data from spliced and unspliced counts
    spliced <- ReadMtx(mtx      = paste0(VelocytoPath,"spliced.mtx.gz"),
                       cells    = paste0(VelocytoPath,"barcodes.tsv.gz"),
                       features = paste0(VelocytoPath, "features.tsv.gz"))

    spliced <- suppressWarnings(CreateSeuratObject(
      counts    = spliced,
      project   = "spliced",
      min.cells = MinCells))

    unspliced <- ReadMtx(mtx      = paste0(VelocytoPath,"unspliced.mtx.gz"),
                         cells    = paste0(VelocytoPath,"barcodes.tsv.gz"),
                         features = paste0(VelocytoPath, "features.tsv.gz"))

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

    cat("✔  \033[1mVelocyto\033[22m nuclear fraction output added\n")

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
      cat("✔  No red blood cells found\n")

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
    saveRDS(hemo_feature_QC, paste0(OutputPath, "RBCQC/", ProjectName, "_RBCQC", ".rds"))

    # Filtering step after visualization, keep those not marked as RBC
    Seurat <- subset(Seurat, RBC == "non-RBC")

    cat("✔  Red blood cells removed\n")

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
    saveRDS(immune_feature_QC, paste0(OutputPath, "IMCQC/", ProjectName, "_IMCQC", ".rds"))

    # Filtering step after visualization, keep those that are immune cells
    Seurat <- subset(Seurat, IMC == "IMC")

    cat("✔  Immune cells isolated\n")
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
  logo <- readPNG("/home/alicen/Projects/limiric/limiric_transparent.png")
  dim <- readPNG("/home/alicen/Projects/limiric/tSNE_right.png")

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


  # Save the final plot
  ggsave(paste0(OutputPath, "/CellQC/", ProjectName, ".png"), plot = final_plot, width = 12, height = 10, dpi = 300)


  cat("✔  \033[3mlimiric\033[23m damaged cell predictions\n")

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
    edresultsDf <- identify_empty_drops(edDf)

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
    dcresultsDf <- identify_damaged_cells(dcDf)

    # Reset output redirection
    sink(NULL)
    sink(NULL, type = "message")

    # Close the connection
    close(tmp_conn)

    # Add DropletQC output to Seurat object metadata
    Seurat$DropletQC <- dcresultsDf$df$cell_status

    # # Rescue cells (dataset specific)
    # if (Organism == "Hsap") {
    #
    #   Seurat$DropletQC <- ifelse(Seurat$DropletQC == "empty_droplet" & Seurat$nf >= 0.02, "cell", Seurat$DropletQC)
    #   Seurat$DropletQC <- ifelse(Seurat$DropletQC == "empty_droplet" & Seurat$nCount_RNA >= 2000, "cell", Seurat$DropletQC)
    # }


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

      cat(paste0("✔  No agreement detected between ", length_damaged_DropletQC, " \033[1mDropletQC\033[22m and ", length_damaged_limiric, " \033[3mlimiric\033[23m damaged cells\n"))

    }

    else

    {

      seurat_damaged <- subset(Seurat, quality == "discard")
      damaged <- length(Cells(seurat_damaged))
      DamagedPercent <- ( damaged / InitialCells ) * 100

      cat("✔  \033[1mDropletQC\033[22m agreement calculated\n")

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
    ggsave(paste0(OutputPath, "/DropletQC/", ProjectName, "_clusters.png"), plot = cluster_plot, width = 5, height = 3, dpi = 300)


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
    clusters <- readPNG(paste0(OutputPath, "/DropletQC/", ProjectName, "_clusters.png"))

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

    # Save the final plot
    ggsave(paste0(OutputPath, "/DropletQC/", ProjectName, ".png"), plot = final_plot, width = 12, height = 10, dpi = 300)


    # Finalise output
    # Save a list of barcodes with QC annnotations (cell, agreed, DropletQC, limiric)
    cleanCells <- data.frame(barcode = rownames(limiric@meta.data), QCannotation = limiric@meta.data$QC)
    write.csv(cleanCells, file = paste0(OutputPath, "/DropletQC/", ProjectName, "_barcodes.csv"), row.names = FALSE)

    # Rename column
    Seurat$limiric.DropletQC <- Seurat$QC

  }

  #  Clean & save output #####


  # Filter cells for saving the Seurat object itself (but only if specified)

  # Filter based on agreement
  if (FilterOutput == TRUE & DropletQC == TRUE) {

    # Save a list of barcodes with limiric annotations (cell, damaged)
    cleanCells <- data.frame(barcode = rownames(Seurat@meta.data), limiric = Seurat@meta.data$limiric)
    write.csv(cleanCells, file = paste0(OutputPath, "/Filtered/", ProjectName, "_barcodes.csv"), row.names = FALSE)

    Seurat <- subset(Seurat, QC != "agreed")

    # Filter unnecessary columns
    columns <- c("nf", "ptprc.percent", "RBC", "IMC", "quality", "nf", "hbb.percent",  "QC", "DropletQC", "limiric.mi", "limiric.ri", "limiric.complexity", "limiric", "limiric.DropletQC")

    for (column in columns){
      if (column %in% colnames(Seurat@meta.data)) {
        Seurat[[column]] <- NULL }
    }


    # Save the filtered Seurat object
    saveRDS(Seurat, paste0(OutputPath, "Filtered/", ProjectName, "_agreement_filtered.rds"))

  }

  # Filter based only on limiric
  if (FilterOutput == TRUE & DropletQC != TRUE) {

    # Save a list of barcodes with limiric annotations (cell, damaged)
    cleanCells <- data.frame(barcode = rownames(Seurat@meta.data), limiric = Seurat@meta.data$limiric)
    write.csv(cleanCells, file = paste0(OutputPath, "/Filtered/", ProjectName, "_barcodes.csv"), row.names = FALSE)

    # Filter damaged cells only according to limiric estimations
    Seurat <- subset(Seurat, limiric == "cell")

    # Filter unnecessary columns
    columns <- c("nf", "ptprc.percent", "RBC", "IMC", "quality", "nf", "hbb.percent",  "QC", "DropletQC", "limiric.mi", "limiric.ri", "limiric.complexity", "limiric", "limiric.DropletQC")

    for (column in columns){
      if (column %in% colnames(Seurat@meta.data)) {
        Seurat[[column]] <- NULL }
    }


    # Save the filtered Seurat object
    saveRDS(Seurat, paste0(OutputPath, "Filtered/", ProjectName, "_filtered.rds"))

  }

  if (FilterOutput == FALSE) {

    # Save object without removing meta data columns (user can visualise for themsleves)
    saveRDS(Seurat, paste0(OutputPath, "Filtered/", ProjectName, "_unfiltered.rds"))

  }

  cat("✔  \033[3mlimiric\033[23m analysis complete.\n\n")

  return(list(
    SeuratObject = Seurat
  ))

}



# Full function ####
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
    if (is.null(SeuratInput)) { SeuratInput = FALSE}
    if (is.null(MinCells)) { MinCells = 0 }
    if (is.null(SoupX)) { SoupX = FALSE}
    if (is.null(DropletQC)) { DropletQC = FALSE}
    if (is.null(FilterRBC)) { FilterRBC = TRUE}
    if (is.null(IsolateCD45)) { IsolateCD45 = FALSE}
    if (is.null(FilterOutput)) { FilterOutput = TRUE}
    if (is.null(Organism)) { Organism = "Hsap"}

    # Run single sample (default) using the core function
    limiricCore(

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

    # Red blood cell QC
    RBCQC_path <- paste0(OutputPath, "RBCQC/", ProjectName, "_RBCQC", ".rds")
    RBCQC      <- readRDS(RBCQC_path)
    ggsave(paste0(OutputPath, "RBCQC/", ProjectName, "_RBCQC", ".png"), plot = RBCQC, width = 5, height = 5, dpi = 300)

    file.remove(RBCQC_path)

    if (IsolateCD45 == TRUE) {

      # And immune QC
      IMCQC_path <- paste0(OutputPath, "IMCQC/", ProjectName, "_IMCQC", ".rds")
      IMCQC      <- readRDS(IMCQC_path)
      ggsave(paste0(OutputPath, "IMCQC/", ProjectName, "_IMCQC", ".png"), plot = IMCQC, width = 5, height = 5, dpi = 300)
      file.remove(IMCQC_path)

    }


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
      if (is.null(SeuratInput)) { SeuratInput = FALSE}
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

        results[[ProjectName]] <- temp_result$SeuratObject

      },

      # If error in one sample, print message and continue with next sample
      error = function(e) {
        cat(paste(e$message, "\n", "Error in processing", ProjectName, "\n", "Likely low cell number\n\n"))

      })

    }

    # Collate output QC plots ####
    # RBC
    rds_dir   <- paste0(OutputPath, "RBCQC/") # Define the directory containing the RDS files
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

    for (i in seq(1, num_plots, by = plots_per_page)) {
      end_index <- min(i + plots_per_page - 1, num_plots)
      current_plots <- plots[i:end_index]

      # Calculate the number of rows needed for the current set of plots
      num_current_plots <- length(current_plots)
      nrow <- ceiling(num_current_plots / 5)

      # Save the current set of plots as a PNG
      png_path <- paste0(OutputPath, "RBCQC/", ProjectName, "_RBCQC_", page_num, ".png")
      save_plot_grid(current_plots, png_path, ncol = 5, nrow = nrow)

      page_num <- page_num + 1

    }

    # Delete all rds files (only want output image)
    file.remove(rds_files)


    if (!is.null(IsolateCD45) & IsolateCD45 == TRUE) {

      # Collate output IMC QC plots ####
      # IMC
      rds_dir   <- paste0(OutputPath, "IMCQC/") # Define the directory containing the RDS files
      rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
      plots     <- lapply(rds_files, readRDS) # Read all RDS files into a list of ggplot objects

      # Loop through the plots & calculate # plot pages
      num_plots <- length(plots)
      plots_per_page <- 15
      page_num <- 1

      for (i in seq(1, num_plots, by = plots_per_page)) {
        end_index <- min(i + plots_per_page - 1, num_plots)
        current_plots <- plots[i:end_index]

        # Calculate the number of rows needed for the current set of plots
        num_current_plots <- length(current_plots)
        nrow <- ceiling(num_current_plots / 5)

        # Save the current set of plots as a PNG
        png_path <- paste0(OutputPath, "IMCQC/", ProjectName, "_IMCQC_", page_num, ".png")
        save_plot_grid(current_plots, png_path, ncol = 5, nrow = nrow)

        page_num <- page_num + 1
      }

      # Delete all rds files (only want output image)
      file.remove(rds_files)
    }

    return(results)

  }

}
