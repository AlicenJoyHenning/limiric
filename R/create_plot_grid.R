#' Save a Grid of Plots
#'
#' @name create_plot_grid
#'
#' @description Helper function to save a grid of plots with a specified number of columns and rows.
#'
#' @param plots A list of 'ggplot2' objects to be arranged in a grid.
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
#' 
#' @examples
#' \dontrun{
#' # Example usage:
#' # Assuming `plot_list` is a list of ggplot objects.
#' plot_list <- list(
#'   ggplot(mtcars, aes(x = mpg, y = wt)) + geom_point(),
#'   ggplot(mtcars, aes(x = hp, y = qsec)) + geom_point(),
#'   ggplot(mtcars, aes(x = drat, y = wt)) + geom_point()
#' )
#' file_path <- "path/to/save/plot_grid.png"
#' nrow <- 2
#' ncol <- 2
#' create_plot_grid(plot_list, file_path, nrow, ncol)
#' }

create_plot_grid <- function(plots, file_path, nrow, ncol = 5) {
  plot_grid <- plot_grid(plotlist = plots, ncol = ncol, nrow = nrow)
  ggsave(file_path, plot = plot_grid, width = ncol * 5, height = nrow * 5, dpi = 300)
}
