#' Identify empty droplets
#'
#' @name identify_empty_droplets
#'
#' @description powellgenomicslab/DropletQC function is used to identify a suitable nuclear fraction
#'   cut-off point to guide the identification of empty droplets. To do this it
#'   calculates the kernel density estimate of the input nuclear fraction scores
#'   and identifies the trough after the first peak, which is assumed to
#'   represent the population of empty droplets.
#'
#' @param nf_umi data frame, containing two columns; the nuclear fraction
#'   estimates in the first column and the total UMI count for each barcode in
#'   the second column
#' @param nf_rescue numeric, a rescue parameter defining a minimum nuclear
#'   fraction score between zero and one. This is used in combination with
#'   `umi_rescue` to identify cells that were misidentified as empty droplets
#' @param umi_rescue integer, a rescue parameter defining a minimum UMI count.
#'   This is used in combination with `nf_rescue` to identify cells that were
#'   misidentified as empty droplets
#' @param plot_name character, if provided a plot will be saved with the
#'   provided name
#' @param plot_path character, if provided a plot will be saved to the provided
#'   path
#' @param plot_width numeric, plot width in cm
#' @param plot_height numeric, plot height in cm
#' @param pdf_png character, either "png" or "pdf"
#'
#' @return data frame, the original data frame is returned plus an additional
#'   column identifying each barcode as a "cell" or "empty_droplet"
#' @export
#'
#' @examples
#' \dontrun{
#' data("qc_examples")
#' gbm <- qc_examples[qc_examples$sample=="GBM",]
#' gbm.ed <- gbm[,c("nuclear_fraction_droplet_qc","umi_count")]
#' gbm.ed <- identify_empty_drops(nf_umi = gbm.ed)
#' head(gbm.ed)
#' table(gbm.ed$cell_status)
#' }

identify_empty_droplets <-
  function(nf_umi,
           nf_rescue = 0.05,
           umi_rescue = 1000,
           plot_name = NULL,
           plot_path = NULL,
           plot_width=18,
           plot_height=13,
           pdf_png = "png"
  ) {

    ## Check and parse arguments
    if (any(class(nf_umi) == "data.frame")) {

      # Assume nuclear fraction is in the first column
      nf <- unlist(nf_umi[, 1], use.names = FALSE)
      # Assume UMI counts are in the second column
      umi <- unlist(nf_umi[, 2], use.names = FALSE)

      # Check values are reasonable
      if(any(c(max(nf)>1, min(nf)<0))){
        warning(paste0("The nuclear fraction values provided in the first column of 'nf_umi' should be between 0 and 1, but values outside this range were identified : minimum = ",min(nf),", maximum = ",max(nf)), call.=FALSE)
      }

      if(!all(umi == floor(umi))){
        non_integer_examples <- which(umi != floor(umi))
        if(length(non_integer_examples)>5){
          non_integer_examples <- non_integer_examples[1:5]
        }
        non_integer_examples <- paste(umi[non_integer_examples], collapse = ",")
        warning(paste0("Non-integer values detected in the second column of 'nf_umi' (e.g. ",non_integer_examples,") where umi counts were expected"), call.=FALSE)
      }

      if(max(umi)<100){
        umi
        warning(paste0("The total umi counts provided in the second column of 'nf_umi' appear to be quite low (max = ",max(umi),"), are these the total UMI counts per cell?"), call.=FALSE)
      }

    } else {
      stop(paste0("A data frame should be supplied to the nf_umi argument, but an object of class ",paste(class(nf_umi), collapse = "/")," was provided"), call.=FALSE)
    }

    # Density estimation (automatically chosen bandwidth)
    kdde_0 <- ks::kdde(x = nf, deriv.order = 0)
    kdde_0 <- data.frame(estimate = kdde_0[["estimate"]],
                         eval.points = kdde_0[["eval.points"]])

    # Density derivative estimation (automatically chosen bandwidth, but different
    # from kdde_0!)
    kdde_1 <- ks::kdde(x = nf, deriv.order = 1)
    kdde_1 <- data.frame(estimate = kdde_1[["estimate"]],
                         eval.points = kdde_1[["eval.points"]])

    # Find point to place cut-off between empty droplets and cells
    gradient_sign <- rle(kdde_1[["estimate"]]>0)
    nf_cutoff <- kdde_1[["eval.points"]][sum(gradient_sign[["lengths"]][1:2])]

    ## Check if there is more than one peak
    if(length(gradient_sign$values)<4){
      warning(paste0("Could not detect more than one peak in the nuclear fraction distribution. There may not be any empty droplets present. We suggest visualising the density estimation (include_plot=TRUE)."), call.=FALSE)
    }

    # Label cells
    nf_umi$cell_status <- "cell"
    nf_umi$cell_status[nf < nf_cutoff] <- "empty_droplet"

    # Rescue miscalled cells
    nf_umi$cell_status[nf > nf_rescue &  umi > umi_rescue] <- "cell"

    return(nf_umi)
  }
