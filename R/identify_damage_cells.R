# This file includes code derived from the DropletQC package.
#
# Copyright (c) 2021 Walter Muskovic
#
# The code is licensed under the MIT license. See the LICENSE file for more details.

#' Identify damage cells
#'
#' @name identify_damage_cells
#'
#' @description This function uses a combination of the cell UMI counts and the
#'   nuclear fraction score to assign each cell one of two values; "cell" or
#'   "damaged_cell". This is based on the idea that damaged cells have a lower
#'   UMI count and higher nuclear fraction than whole cells. The expected input
#'   is a data frame with four columns. The first three columns should contain;
#'   the nuclear fraction score, total UMIs and a character vector describing
#'   each cell as "cell" or "empty_droplet". This is the format output by the
#'   `identify_empty_drops` function. The fourth column should be a character
#'   vector with user-assigned cell types. Internally, the provided data frame
#'   is split by cell type and a Gaussian mixture model with a maximum of two
#'   components is fit to the umi counts and nuclear fraction scores. The
#'   parameters of the model are estimated using expectation maximisation (EM)
#'   with the `mclust` package. The best model is selected using the Bayesian
#'   Information Criterion (BIC). The two populations (cells and damaged cells)
#'   are assumed to have equal variance (mclust model name "EEI").
#'   EM Model: It’s used to find and separate multiple groups or clusters within
#'   the data. It’s not just about finding points that don’t fit in; it’s about
#'   identifying distinct groups that may exist in the data.
#'
#' @param nf_umi_ed_ct data frame, with four columns. The first three columns
#'   should match the output from the `identify_empty_drops` function. The
#'   fourth column should contain cell type names.
#' @param nf_sep numeric, the minimum separation of the nuclear fraction score
#'   required between the cell and damage cell populations
#' @param umi_sep_perc numeric, this is the minimum percentage of UMIs which the
#'   damaged cell population is required to have compared to the cell
#'   population. For example, if the mean UMI of the distribution fit to the
#'   whole cell population is 10,000 UMIs, the mean of the distribution fit to
#'   the damaged cell population must be at less than 7,000 UMIs if the umi_sep
#'   parameter is 30 (%)
#' @param verbose logical, whether to print updates and progress while fitting
#'   with EM
#'
#' @return list, of length two. The first element in the list contains a data
#'   frame with the same dimensions input to the `nf_umi_ed_ct` argument, with
#'   "damaged_cell" now recorded in the third column.
#'
#' @export
#'
#' @importFrom ks kde hpi
#' @importFrom mclust mclustBIC Mclust
#' @import stats
#'
#' @keywords internal

utils::globalVariables(c("nf_umi"))

identify_damage_cells <- function(nf_umi_ed_ct,
                                  nf_sep       = 0.15,  # Nuclear fraction separation threshold
                                  umi_sep_perc = 50,    # UMI counts percentage less than cell
                                  verbose      = TRUE   # Print progress messages
) {
  # Check if 'nf_umi_ed_ct' is a data frame
  if (any(class(nf_umi_ed_ct) == "data.frame")) {

    # Ensure the data frame has exactly four columns
    if (ncol(nf_umi_ed_ct) != 4) {
      stop("nf_umi_ed_ct should be a data frame with four columns", call. = FALSE)
    }


    nf  <- unlist(nf_umi_ed_ct[, 1], use.names = FALSE) # Extract the nuclear fraction (1st column)
    umi <- unlist(nf_umi_ed_ct[, 2], use.names = FALSE) # Extract UMI counts (2nd column)
    ed  <- unlist(nf_umi_ed_ct[, 3], use.names = FALSE) # Extract cell status ("cell" or "empty_droplet", 3rd column)
    ct  <- unlist(nf_umi_ed_ct[, 4], use.names = FALSE) # Extract cell types (4th column)


    # Validate nuclear fraction values are between 0 and 1
    if (any(c(max(nf) > 1, min(nf) < 0))) {
      warning(paste0("Nuclear fraction values should be between 0 and 1, found range: ", min(nf), " to ", max(nf)), call. = FALSE)
    }

    # Ensure UMI counts are integers
    if (!all(umi == floor(umi))) {
      non_integer_examples <- which(umi != floor(umi))
      if (length(non_integer_examples) > 5) {
        non_integer_examples <- non_integer_examples[1:5]
      }
      non_integer_examples <- paste(umi[non_integer_examples], collapse = ",")
      warning(paste0("Non-integer UMI counts detected (e.g., ", non_integer_examples, ")"), call. = FALSE)
    }

    # Warn if UMI counts appear too low
    if (max(umi) < 100) {
      warning(paste0("UMI counts appear low (max = ", max(umi), "), ensure these are total UMI counts per cell"), call. = FALSE)
    }

    # Validate cell status values
    if (!all(unique(ed) %in% c("cell", "empty_droplet"))) {
      ed_output <- unique(ed)
      if (length(ed_output) > 5) {
        ed_output <- ed_output[1:5]
      }
      ed_output <- paste(ed_output, collapse = ",")
      warning(paste0("Unexpected values in the third column: ", ed_output), call. = FALSE)
    }

    # Print unique cell types if verbose is TRUE
    if (verbose) {
      ct <- unique(ct)
      ct <- paste(ct, collapse = ",")
      print(paste0("Provided cell types: ", ct))
    }

  } else {
    stop(paste0("Expected a data frame for 'nf_umi_ed_ct', but received: ", paste(class(nf_umi_ed_ct), collapse = "/")), call. = FALSE)
  }

  # Prepare data for EM (Expectation-Maximization) model fitting
  em.data <- data.frame(
    nf  = unlist(nf_umi_ed_ct[, 1], use.names = FALSE),        # Nuclear fraction
    umi = log10(unlist(nf_umi_ed_ct[, 2], use.names = FALSE)), # Log-transformed UMI counts
    ct  = unlist(nf_umi_ed_ct[, 4], use.names = FALSE)         # Cell types
  )
  row.names(em.data) <- 1:nrow(em.data)

  # Filter out empty droplets
  em.data <- em.data[nf_umi_ed_ct[, 3] == "cell", ]

  # Split data by cell type
  em.data.ct <- split(em.data, em.data$ct)

  # Fit EM models (assess_EM) to each cell type
  if (verbose) {
    print("Fitting models with EM")
  }
  em_mods <- lapply(em.data.ct, function(x) mclust::Mclust(data = x[, 1:2], G = 1:2, modelNames = "EEI", verbose = verbose))

  # Assess EM models and assign barcodes as "cell" or "damaged_cell"
  em_mods_assessed <- lapply(em_mods, assess_EM, nf_thresh = nf_sep, umi_thresh = umi_sep_perc)

  # Update the original data frame with damaged cell info
  names(em_mods_assessed) <- NULL
  em_mods_assessed <- do.call(rbind, em_mods_assessed)
  nf_umi_ed_ct$cell_status[as.integer(row.names(em_mods_assessed))] <- em_mods_assessed$classification

  # Return results as a list containing the updated data frame
  return(list(df = nf_umi_ed_ct, plots = NULL))
}

