#' Identify damaged cells with DamageDetective
#'
#' @name DamageDetective_calculation
#'
#' @description This helper function calculates the 'DamageDetective' metrics for a given 'Seurat' object and organism.
#' It identifies mitochondrial and ribosomal genes, reduces the 'Seurat' object based on these genes,
#' and calculates various metrics including mitochondrial and ribosomal expression complexity.
#'
#' @param organism A string representing the organism ("Hsap", "Mmus").
#' @param Seurat A 'Seurat' object containing the single-cell RNA-seq data.
#' @param resolution Numeric between 0 and 1.6 describing cluster division. Default is 1.
#' @param cluster_ranks Numeric describing the number of top ranking clusters to be included as damaged cells. Default 1.
#' @param annotations A data frame containing gene annotations.
#' @param initial_cells An integer representing the initial number of cells.
#' @param project_name A string representing the name of the project, used for plot titles.
#' @param output_path A string representing the path to save the output plots.
#'
#' @return A list containing the updated 'Seurat' object and the reduced 'Seurat' object.
#'
#' @import cowplot
#' @importFrom dplyr %>% pull group_by summarise mutate arrange slice
#' @import ggplot2
#' @importFrom png readPNG
#' @import Seurat
#' @importFrom utils globalVariables
#' @importFrom stats quantile
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # Load test data
#' data("test_data", package = "DamageDetective")
#' data("human_annotations", package = "DamageDetective")
#'
#' result <- DamageDetective_calculation(organism = "Hsap",
#'                               Seurat = test_data,
#'                               annotations = human_annotations,
#'                               resolution = 0.1,
#'                               cluster_ranks = 1,
#'                               initial_cells = 1000,
#'                               project_name = "Test",
#'                               output_path = tempdir())
#' print(result)
#' }

utils::globalVariables(c("mt_plot", "complexity_plot", "rb_plot",
                         "cluster_plot", "logo", "DamageDetective", "Cells"))

DamageDetective_calculation <- function(organism,
                                Seurat,
                                resolution,
                                cluster_ranks,
                                annotations,
                                initial_cells,
                                project_name,
                                output_path
) {
  # Define mitochondrial and ribosomal low dimensional gene space ------------------------------------

  if (organism == "Hsap") {

    # Get gene annotations for mitochondrial (MT) & ribosomal (RB) genes
    mt_gene_annotations <- annotations[grep("^MT-", annotations$gene_name, perl=TRUE),]
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

  if (organism == "Mmus") {

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


  # Dimensionality reduction step ------------------------------------

  # Reduce based on mt & rb genes only (this occurs in a separate Seurat object (DamageDetective))
  DamageDetective <- subset(Seurat, features = intersect(mt_rb_genes, rownames(Seurat@assays$RNA)))

  DamageDetective <- NormalizeData(DamageDetective, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE) %>%
    RunTSNE(dims = 1:10, verbose = FALSE, check_duplicates = FALSE)

  # Calculate mitochondrial & ribosomal QC metrics ------------------------------------

  # Annotations to the reduced object based on the unreduced seurat object (all genes)
  DefaultAssay(DamageDetective) <- "RNA"
  DefaultAssay(Seurat)  <- "RNA"

  # Define mitochondrial expression
  DamageDetective$mt.percent <- PercentageFeatureSet(
    object   = Seurat,
    features = intersect(mt_genes, rownames(Seurat@assays$RNA)),
    assay    = "RNA"
  )

  # Transfer to actual Seurat object
  Seurat$DamageDetective.mi <- DamageDetective$mt.percent

  # Define ribosomal expression
  DamageDetective$rb.percent <- PercentageFeatureSet(
    object   = Seurat,
    features = intersect(rb_genes, rownames(Seurat@assays$RNA)),
    assay    = "RNA"
  )

  # Transfer to actual Seurat object
  Seurat$DamageDetective.ri <- DamageDetective$rb.percent

  # Identify damaged cell cluster using QC metric averages for each cluster ------------------------------------

  # Automatically find the damaged cell population
  best_cluster <- DamageDetective@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      avg_mt_percent = mean(mt.percent, na.rm = TRUE),
      avg_rb_percent = mean(rb.percent, na.rm = TRUE)
    ) %>%
    mutate(
      rank_mt = rank(-avg_mt_percent),
      rank_rb = rank(avg_rb_percent),
      combined_rank = rank_mt + rank_rb
    ) %>%
    arrange(combined_rank) %>%
    slice(1:cluster_ranks) %>%
    pull(seurat_clusters) %>%
    as.character()

  # Label all cells belonging to this cluster as "damaged"
  DamageDetective$seurat_clusters <- ifelse(DamageDetective$seurat_clusters %in% best_cluster, 'damaged', DamageDetective$seurat_clusters)


  # Calculate complexity score for each cell ----

  # Min-max scaling function
  min_max_scale <- function(x) {
    return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }

  # Apply min-max scaling for mt & rb percent expression and library size
  # Value is out of 3 where closer to 0 is less complex, more likely to be damaged
  DamageDetective$complexity <- ((DamageDetective$rb.percent / 100 ) +
                                (1 - (DamageDetective$mt.percent / 100 )) + # Inverse for mt.percent, highe it is the less complex cells are
                                 min_max_scale(DamageDetective$nFeature_RNA)) / 3


  # Find the median complexity score for all undamaged cells
  clusters <- table(DamageDetective$seurat_clusters)

  if (length(clusters) >= 2){
  undamaged_cells <- subset(DamageDetective, seurat_clusters != "damaged")
  complexity_threshold <- quantile(undamaged_cells$complexity, probs = 0.05, na.rm = TRUE)

  } else {
  complexity_threshold <- quantile(DamageDetective$complexity, probs = 0.05, na.rm = TRUE)
  }

  # Use this value as minimum threshold for retaining damaged label
  DamageDetective$seurat_clusters <- ifelse(

     # First, looking at cells in damage cluster
     DamageDetective$seurat_clusters == "damaged" &

     # Next, retain if they have sufficiently low complexity
     DamageDetective$complexity <= complexity_threshold,

     "damaged", 0)


  # Add Cell QC meta data to object
  DamageDetective$MtRb <- ifelse(DamageDetective$seurat_clusters == "damaged", "damaged", "cell")

  # Calculated percentage of damaged cells
  damaged <- subset(DamageDetective, MtRb == "damaged")
  damagedCells <- length(Cells(damaged))
  damaged_percent <- (damagedCells / initial_cells) * 100

  # Add cluster #s to actual seurat object
  Seurat$DamageDetective <- DamageDetective$MtRb


  # Visualise  clusters coloured by QC metrics ------------------------------------

  mt_plot <- FeaturePlot(DamageDetective, features = c("mt.percent"), cols = c("#E1E1E1", "#0073CF"), pt.size = 1) +
    NoAxes() + labs(caption = "Mitochondrial gene expression") +
    theme(
      plot.title   = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))

  rb_plot <- FeaturePlot(DamageDetective, features = c("rb.percent"), cols = c("#E1E1E1", "#0073CF"), pt.size = 1) +
    NoAxes() + labs(caption = "Ribosomal gene expression") +
    theme(
      plot.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))

  cluster_plot <- DimPlot(
    DamageDetective, pt.size = 1, group.by = "MtRb", cols = c("cell" = "grey", "damaged" = "maroon")) +
    labs(caption = expression("Damaged cells identified by " * italic("DamageDetective"))) + NoAxes() +
    theme(
      plot.title = element_blank(),
      plot.caption = element_text(hjust = 0.5, size = 16),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1)
    )

  complexity_plot <- FeaturePlot(DamageDetective, features = c("nFeature_RNA"), cols = c("#E1E1E1", "#0073CF"), pt.size = 1) +
    NoAxes() + labs(caption = "Library size") +
    theme(
      plot.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))

  title <- project_name
  label <- paste("Estimated damaged cells: ", round(damaged_percent, 2), "%, ", initial_cells, " cells")
  dim <- readPNG(system.file("extdata", "tSNE.png", package = "DamageDetective"))


  DamageDetective_plot <- (mt_plot | complexity_plot) / (rb_plot | cluster_plot)

  title <- ggdraw() + draw_label(project_name, fontface = 'bold', x = 0.45, y = 0.1, hjust = 0.5, size = 20)
  subtitle <- ggdraw() + draw_label(paste("Estimated", round(damaged_percent, 2), "% damaged of ", initial_cells, " cells"), x = 0.45, hjust = 0.5, size = 16)

  dim_plot <- ggdraw() + draw_image(dim, x = 0.54, y = 0.5, width = 0.8, height = 0.8)

  final_plot <- plot_grid(
    title, subtitle, DamageDetective_plot, dim_plot,
    ncol = 1,
    rel_heights = c(0.1, 0.1, 1, 0.2)
  )

  # Set the background to white for the entire plot
  final_plot <- final_plot + theme(plot.background = element_rect(fill = "white", color = "white"))

  # Save the final plot
  ggsave(file.path(output_path, "/CellQC/", paste0(project_name, ".png")), plot = final_plot, width = 12, height = 10, dpi = 300)


  return(list(Seurat = Seurat, DamageDetective = DamageDetective))

}
