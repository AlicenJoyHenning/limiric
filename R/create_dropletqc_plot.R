#' Create DropletQC plot output
#'
#' @name create_dropletqc_plot
#'
#' @description This helper function generates a cluster and scatter plot
#' for viewing labels in reduced space and scatter plot of nf vs UMI_count.
#' The result is saved as a PNG.
#'
#' @param Seurat A Seurat object containing the single-cell RNA-seq data.
#' @param limiric A Seurat object containing the limiric data.
#' @param output_path A string representing the path to save the plots.
#' @param project_name A string representing the name of the project, used as the plot title.
#' @param damaged_percent A numeric value representing the percentage of damaged cells.
#' @param initial_cells An integer representing the initial number of cells.
#'
#' @return A list containing the Seurat object and the limiric object.
#'
#' @import cowplot
#' @importFrom dplyr %>% pull group_by summarise mutate arrange slice case_when
#' @import ggplot2
#' @importFrom png readPNG
#' @import Seurat
#' @importFrom utils write.csv
#'
#' @export
#'
#' @keywords internal

utils::globalVariables(c("QC", "nf", "nCount_RNA"))

create_dropletqc_plot <- function(Seurat,
                                  limiric,
                                  output_path,
                                  project_name,
                                  damaged_percent,
                                  initial_cells
) {
  # View labels in reduced space ------------------------------------

  cluster_plot <- DimPlot(limiric, pt.size = 1,
                          group.by = "QC",
                          cols = c("cell"       = "#D5D5D5",
                                   "agreed"     = "#506BB0",
                                   "droplet_qc" = "#D98C25",
                                   "limiric"    = "#A4BAF2")) + NoAxes() +
    theme(plot.title   = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

  # Save the cluster plot
  ggsave(file.path(output_path, "droplet_qc", paste0(project_name, "_clusters.png")),
         plot = cluster_plot, width = 5, height = 3, dpi = 300)



  # View labels in scatter plot nf vs UMI_count (droplet_qc metric) ------------------------------------

  df <- Seurat@meta.data
  df$nf <- as.numeric(df$nf)
  df$nCount_RNA <- as.numeric(df$nCount_RNA)

  scatter <- ggplot(df, aes(x = nf, y = nCount_RNA, color = QC)) +
    geom_point(size = 1.2) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "agreed" = "#506BB0", "droplet_qc" = "#D98C25", "limiric" = "#A4BAF2")) +
    coord_trans(y = "log10") +
    xlab("nf") + ylab("log10(UMI)") +
    labs(caption = expression("Agreement between " * italic("droplet_qc") * " and " * italic("limiric") * " estimations")) +
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
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_line(colour = "black")
    )

  # FANCY PLOT COMBINING
  clusters <- readPNG(file.path(output_path, "droplet_qc", paste0(project_name, "_clusters.png")))

  # Create the title and subtitle
  title <- ggdraw() + draw_label(project_name, fontface = 'bold', x = 0.5, y = 0.1, hjust = 0.5, size = 20)
  subtitle <- ggdraw() + draw_label(paste("Agreed estimated", round(damaged_percent, 2), "% damaged of ", initial_cells, " cells"), x = 0.5, hjust = 0.5, size = 16)

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
  ggsave(file.path(output_path, "droplet_qc", paste0(project_name, ".png")), plot = final_plot, width = 12, height = 10, dpi = 300)

  invisible()

}
