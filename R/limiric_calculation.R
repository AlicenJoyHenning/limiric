#' Identify damaged cells with limiric
#'
#' @name limiric_calculation
#'
#' @description This helper function calculates the LIMIRIC metrics for a given Seurat object and organism.
#' It identifies mitochondrial and ribosomal genes, reduces the Seurat object based on these genes,
#' and calculates various metrics including mitochondrial and ribosomal expression, and complexity.
#'
#' @param organism A string representing the organism ("Hsap", "Mmus").
#' @param Seurat A Seurat object containing the single-cell RNA-seq data.
#' @param annotations A data frame containing gene annotations.
#' @param initial_cells An integer representing the initial number of cells.
#' @param project_name A string representing the name of the project, used for plot titles.
#' @param output_path A string representing the path to save the output plots.
#'
#' @return A list containing the updated Seurat object and the reduced Seurat object (limiric).
#'
#' @import cowplot
#' @importFrom dplyr %>% pull group_by summarise mutate arrange slice
#' @import ggplot2
#' @importFrom png readPNG
#' @import Seurat
#' @importFrom utils globalVariables
#'
#' @export
#'
#' @keywords internal

utils::globalVariables(c("mt_plot", "complexity_plot", "rb_plot",
                         "cluster_plot", "logo", "limiric", "Cells"))

limiric_calculation <- function(organism,
                                Seurat,
                                annotations,
                                initial_cells,
                                project_name,
                                output_path
) {
  # Define mitochondrial and ribosomal low dimensional gene space ------------------------------------

  if (organism == "Hsap") {

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

  # Reduce based on mt & rb genes only (this occurs in a separate Seurat object (limiric))
  limiric <- subset(Seurat, features = intersect(mt_rb_genes, rownames(Seurat@assays$RNA)))

  limiric <- NormalizeData(limiric, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE) %>%
    RunTSNE(dims = 1:10, verbose = FALSE, check_duplicates = FALSE)

  # Calculate mitochondrial & ribosomal QC metrics ------------------------------------

  # Annotations to the reduced object based on the unreduced seurat object (all genes)
  DefaultAssay(limiric) <- "RNA"
  DefaultAssay(Seurat)  <- "RNA"

  # Define mitochondrial expression
  limiric$mt.percent <- PercentageFeatureSet(
    object   = Seurat,
    features = intersect(mt_genes, rownames(Seurat@assays$RNA)),
    assay    = "RNA"
  )

  # Transfer to actual Seurat object
  Seurat$limiric.mi <- limiric$mt.percent

  # Define ribosomal expression
  limiric$rb.percent <- PercentageFeatureSet(
    object   = Seurat,
    features = intersect(rb_genes, rownames(Seurat@assays$RNA)),
    assay    = "RNA"
  )

  # Transfer to actual Seurat object
  Seurat$limiric.ri <- limiric$rb.percent

  # Calculate complexity
  matrix <- Seurat::GetAssayData(limiric, layer = "data")
  cat("check")
  complexity_metric <- colSums(matrix > 0)
  limiric$complexity <- complexity_metric

  # Transfer to actual Seurat object
  Seurat$limiric.complexity <- limiric$complexity



  # Identify damaged cells using QC metric averages for each cluster ------------------------------------

  # Automatically find the damaged cell population
  best_cluster <- limiric@meta.data %>%
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
    slice(1) %>%
    pull(seurat_clusters) %>%
    as.character()

  # Label all cells belonging to this cluster as "damaged"
  limiric$seurat_clusters <- ifelse(limiric$seurat_clusters == best_cluster, 'damaged', limiric$seurat_clusters)

  # Add Cell QC meta data to object
  limiric$MtRb <- ifelse(limiric$seurat_clusters == "damaged", "damaged", "cell")

  # Calculated percentage of damaged cells
  damaged <- subset(limiric, MtRb == "damaged")
  damagedCells <- length(Cells(damaged))
  damaged_percent <- (damagedCells / initial_cells) * 100

  # Add cluster #s to actual seurat object
  Seurat$limiric <- limiric$MtRb


  # Visualise  clusters coloured by QC metrics ------------------------------------

  mt_plot <- FeaturePlot(limiric, features = c("mt.percent"), cols = c("#E1E1E1", "#5372B4"), pt.size = 1) +
    NoAxes() + labs(caption = "Mitochondrial gene expression") +
    theme(
      plot.title   = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))

  rb_plot <- FeaturePlot(limiric, features = c("rb.percent"), cols = c("#E1E1E1", "#5372B4"), pt.size = 1) +
    NoAxes() + labs(caption = "Ribosomal gene expression") +
    theme(
      plot.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))

  cluster_plot <- DimPlot(
    limiric, pt.size = 1, group.by = "MtRb", cols = c("cell" = "grey", "damaged" = "#5372B4")) +
    labs(caption = expression("Damaged cells identified by " * italic("limiric"))) + NoAxes() +
    theme(
      plot.title = element_blank(),
      plot.caption = element_text(hjust = 0.5, size = 16),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1)
    )

  complexity_plot <- FeaturePlot(limiric, features = c("complexity"), cols = c("#E1E1E1", "#5372B4"), pt.size = 1) +
    NoAxes() + labs(caption = "Complexity score") +
    theme(
      plot.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
      plot.caption = element_text(hjust = 0.5, size = 16))

  title <- project_name
  label <- paste("Estimated damaged cells: ", round(damaged_percent, 2), "%, ", initial_cells, " cells")
  dim <- readPNG(system.file("extdata", "tSNE.png", package = "limiric"))

  limiric_plot <- plot_grid(mt_plot, complexity_plot, rb_plot, cluster_plot, ncol = 2)
  limiric_plot <- (mt_plot | complexity_plot) / (rb_plot | cluster_plot)

  title <- ggdraw() + draw_label(project_name, fontface = 'bold', x = 0.45, y = 0.1, hjust = 0.5, size = 20)
  subtitle <- ggdraw() + draw_label(paste("Estimated", round(damaged_percent, 2), "% damaged of ", initial_cells, " cells"), x = 0.45, hjust = 0.5, size = 16)

  dim_plot <- ggdraw() + draw_image(dim, x = 0.54, y = 0.5, width = 0.8, height = 0.8)

  final_plot <- plot_grid(
    title, subtitle, limiric_plot, dim_plot,
    ncol = 1,
    rel_heights = c(0.1, 0.1, 1, 0.2)
  )

  # Set the background to white for the entire plot
  final_plot <- final_plot + theme(plot.background = element_rect(fill = "white", color = "white"))

  # Save the final plot
  ggsave(file.path(output_path, "/CellQC/", paste0(project_name, ".png")), plot = final_plot, width = 12, height = 10, dpi = 300)


  return(list(Seurat = Seurat, limiric = limiric))

}
