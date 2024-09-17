#' Create RBC Plot
#'
#' @name create_rbc_plot
#'
#' @description This helper function generates a scatter plot showing the relationship
#' between haemoglobin expression and CD45 expression, with points colored based
#' on whether they are RBCs or non-RBCs.
#'
#' @param rbc_df A data frame containing the RBC data.
#' @param project_name A string representing the name of the project, used as the plot title.
#' @param RBC_percent A numeric value representing the percentage of RBC contamination.
#' @param initial_cells An integer representing the initial number of cells.
#'
#' @return A 'ggplot2' object representing the RBC plot.
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
#' # Create test dataset
#' rbc_data <- data.frame(
#'   hemo.percent = runif(100, min = 0, max = 100),
#'   ptprc.percent = runif(100, min = 0, max = 100),
#'   RBC = sample(c("RBC", "non-RBC"), 100, replace = TRUE)
#' )
#'
#' create_rbc_plot(rbc_df = rbc_data,
#'                 project_name = "test",
#'                 RBC_percent = 15,
#'                 initial_cells = 1000)
#'
#' }

utils::globalVariables(c("RBC"))

create_rbc_plot <- function(rbc_df, project_name, RBC_percent, initial_cells) {

  plot <- ggplot(rbc_df, aes(x = hemo.percent, y = ptprc.percent, color = RBC)) +
    geom_point(size = 0.6) +
    scale_color_manual(values = c("non-RBC" = "grey", "RBC" = "#9AA9FF")) + # "#6765ED")) +
    xlab("Haemoglobin expression") +
    ylab("CD45 expression") +
    labs(title = project_name) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste("Contamination of", round(RBC_percent, 2), "% for", initial_cells, "cells"),
      hjust = 1.1, vjust = 2.4
    ) +
    theme_classic() +
    theme(
      axis.line = element_line(colour = "black"),
      axis.title.x = element_text(hjust = 0.5, vjust = -1, face = "bold", size = 12),
      axis.title.y = element_text(hjust = 0.5, vjust = 3, face = "bold", size = 12),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = -1.2, size = 12),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = -1.2, size = 12),
      legend.title = element_text(face = "bold"),
      legend.position = "below",
      legend.box.background = element_rect(colour = "black"),
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.margin = unit(c(1, 1, 1, 1), 'cm'),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )

  return(plot)
}
