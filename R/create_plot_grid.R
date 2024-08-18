#' Save a Grid of Plots
#'
#' @name create_plot_grid
#'
#' @description Helper function to save a grid of plots with a specified number of columns and rows.
#'
#' @param plots A list of ggplot objects to be arranged in a grid.
#' @param file_path A string representing the file path where the plot grid will be saved.
#' @param ncol An integer representing the number of columns in the plot grid. Default is 5.
#' @param nrow An integer representing the number of rows in the plot grid. Default is 3.
#'
#' @return None. The function saves the plot grid to the specified file path.
#'
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggsave
#'
#' @export
#'
#' @keywords internal

create_plot_grid <- function(plots, file_path, nrow, ncol = 5) {
  plot_grid <- plot_grid(plotlist = plots, ncol = ncol, nrow = nrow)
  ggsave(file_path, plot = plot_grid, width = ncol * 5, height = nrow * 5, dpi = 300)
}
