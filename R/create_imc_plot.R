#' Create Immune Feature QC Plot
#'
#' @name create_imc_plot
#'
#' @description Helper function to generate scatter plot of UMI & gene counts per cell barcode
#' with points colored based on whether they are immune cells (IMC) or non-immune cells.
#'
#' @param imc_df A data frame containing the immune cell data.
#' @param project_name A string representing the name of the project, used as the plot title.
#' @param IMC_percent A numeric value representing the percentage of immune cells.
#'
#' @return A ggplot object representing the immune feature QC plot.
#'
#' @import ggplot2
#' @importFrom grid unit
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # Assuming `imc_data` is a data frame containing the immune cell data.
#' imc_data <- data.frame(
#'   nCount_RNA = runif(100, min = 500, max = 10000),
#'   nFeature_RNA = runif(100, min = 200, max = 3000),
#'   IMC = sample(c("IMC", "non-IMC"), 100, replace = TRUE)
#' )
#' project_name <- "ImmuneFeatureQC"
#' IMC_percent <- 25.0
#' plot <- create_imc_plot(imc_data, project_name, IMC_percent)
#' print(plot)
#' }

utils::globalVariables(c("IMC", "nCount_RNA", "nFeature_RNA"))

create_imc_plot <- function(imc_df, project_name, IMC_percent) {

  plot <- ggplot(imc_df, aes(x = nCount_RNA, y = nFeature_RNA, color = IMC)) +
    geom_point(size = 0.6) +
    scale_color_manual(values = c("IMC" = "grey", "non-IMC" = "#6765ED")) +
    labs(title = project_name) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste("Percentage of immune cells:", round(IMC_percent, 2), "%"),
      hjust = 1.1, vjust = 2.4
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.margin = unit(c(1, 1, 1, 1), 'cm'),
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
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )

  return(plot)
}
