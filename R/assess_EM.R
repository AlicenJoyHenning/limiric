# This file includes code derived from the DropletQC package.
#
# Copyright (c) 2021 Walter Muskovic
#
# The code is licensed under the MIT license. See the LICENSE file for more details.

#' assess_EM
#'
#' Assess Expectation-Maximization model results to assign cells as damaged or not
#'
#' @name assess_EM
#'
#' @description powellgenomicslab/DropletQC helper function called by
#'  `identify_damaged_cells` not intended for more general use.
#'  The function is designed to detect and label damaged cells based on
#'  their UMI counts and nuclear fraction, relying on the EM model
#'  to differentiate between potential cell populations. The thresholds
#'  (umi_thresh and nf_thresh) ensure that only cells with significant
#'  differences in these metrics are considered damaged.
#'
#' @param em Mclust, result of EM on log10(UMI counts) and nf to assess
#' the likelihood of two cell populations (cells and damaged cells) or just cells.
#' @param umi_thresh numeric, percentage by which damaged_cell UMI counts
#' must be below the cell distribution mean to be classified as damaged cells
#' (e.g., if umi_thresh = 30, damaged_cell UMI mean must be < 70% of cell mean,
#'  otherwise all barcodes are returned as cells).
#' @param nf_thresh numeric, minimum nuclear fraction difference required
#' between cell and damaged_cell distributions to classify as damaged cells.
#'
#' @return data frame, with three columns containing the log10(UMI counts),
#'   nuclear fraction score and the assigned cell status; "cell" or
#'   "damaged_cell"
#'
#' @importFrom dplyr %>% pull group_by summarise mutate arrange slice case_when
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @importFrom utils data write.csv
#'
#' @keywords internal
#' 
#' @examplesIf interactive()
#' # Example usage:
#' # Assuming `em_result` is the result of an EM model on log10(UMI counts) and nf
#' # and you have set thresholds for UMI and nuclear fraction.
#' em_result <- Mclust(data = your_data)  # Replace `your_data` with actual data
#' umi_threshold <- 30
#' nf_threshold <- 0.1
#' result <- assess_EM(em_result, umi_threshold, nf_threshold)
#' print(result)
#'
#' @export

utils::globalVariables(c("nf_means", "umi_means", "nf_check",
                         "umi_check", "check_1", "check_2",
                         "check_3", "damaged_cells"))

assess_EM <- function(em, umi_thresh, nf_thresh) {

  # This function assigns each cell barcode as "cell" or "damaged_cell"
  # based on the results of the EM model. The assignment is done through
  # the following sequential procedure:

  # Step 1 ------------------------------------
  # Check if the EM model identified two distinct distributions (G=2).
  # If only one distribution is found, classify all barcodes as "cell".
  check_1 <- em$G == 2

  # Step 2 ------------------------------------
  # If two distributions are identified, check if the distribution with
  # the higher nuclear fraction mean also has a lower UMI mean.
  # If true, this suggests that the population with the lower UMI count
  # and higher nuclear fraction might be damaged cells, and we proceed
  #  to the next step. Otherwise, classify all barcodes as "cell".
  if (check_1) {
    nf_means  <- em$parameters$mean["nf", ]  # Mean nuclear fraction for each distribution
    umi_means <- em$parameters$mean["umi", ] # Mean UMI count for each distribution

    # Check if the distribution with the highest nuclear fraction
    # also has the lowest UMI count.
    check_2   <- umi_means[which.max(nf_means)] < umi_means[which.min(nf_means)]
  }

  else {check_2 <- FALSE}

  # Step 3 ------------------------------------
  # If the previous conditions are met, further validate the classification by applying thresholds.
  #    - Check if the difference in nuclear fraction means between the two
  #      distributions is greater than the nf_thresh.
  #    - Check if the UMI count of the "damaged_cell" population is at least
  #      umi_thresh percent lower than the "cell" population.
  # If both conditions are satisfied, classify the relevant barcodes as "damaged_cell"; otherwise, classify all as "cell".

  if (check_2) {
    # Validate nuclear fraction threshold
    nf_check  <- nf_means[which.max(nf_means)] - nf_means[which.min(nf_means)] > nf_thresh

    # Validate UMI threshold
    damaged_cell_umi <- 10^umi_means[which.max(nf_means)] # UMI mean for "damaged_cell"
    cell_umi  <- 10^umi_means[which.min(nf_means)]        # UMI mean for "cell"
    umi_check <- damaged_cell_umi < (cell_umi - cell_umi * (umi_thresh / 100))

    # Confirm that both thresholds are met
    check_3 <- all(nf_check, umi_check)
  } else {
    check_3 <- FALSE
  }

  # Final classification of barcodes based on the checks
  em_classification <- data.frame(
    nf  = em$data[,"nf"],               # Nuclear fraction for each barcode
    umi = em$data[,"umi"],              # UMI count for each barcode
    classification = em$classification  # Initial classification from EM
  )
  row.names(em_classification) <- row.names(em$data)

  # If all checks are passed, assign the relevant barcodes to "damaged_cell"
  # and the others to "cell". Otherwise, classify all as "cell".
  if (all(check_1, check_2, check_3)) {
    damaged_cells <- em_classification$classification == which.max(em$parameters$mean["nf", ])
    em_classification$classification[damaged_cells] <- "damaged_cell"
    em_classification$classification[!damaged_cells] <- "cell"
  } else {
    em_classification$classification <- "cell"
  }

  return(em_classification)
}

