#' TRRident: T Cell Entry Region Identification in Spatial Transcriptomics
#'
#' A computational framework for identifying T cell entry regions (TRR)
#' in spatial transcriptomics data using Ising models with impurities.
#'
#' @section Main Functions:
#'
#' \strong{Data Preparation:}
#' \itemize{
#'   \item \code{\link{data_preparation_TRR}} - Main data preparation function
#'   \item \code{\link{read.spatial.dat}} - Load 10x Visium data
#'   \item \code{\link{plot_seurat_clusters}} - Visualize Seurat clusters
#' }
#'
#' \strong{Ising Model Estimation:}
#' \itemize{
#'   \item \code{\link{estimate_J1_dat}} - Estimate J1 parameter
#'   \item \code{\link{find_optimal_mu_Q}} - Variational inference for J
#'   \item \code{\link{plot_J2}} - Estimate J2 (impurity coupling)
#' }
#'
#' \strong{Line Management:}
#' \itemize{
#'   \item \code{\link{get_lines}} - Decompose 2D data into 1D lines
#'   \item \code{\link{plot_lines}} - Visualize lines
#'   \item \code{\link{get_cps_and_save}} - Detect change points
#'   \item \code{\link{calculate_slopes}} - Calculate segment slopes
#' }
#'
#' \strong{TRR Identification:}
#' \itemize{
#'   \item \code{\link{gibbs_sampling_new}} - Gibbs sampling
#'   \item \code{\link{calculate_final_prob}} - Calculate posterior probabilities
#'   \item \code{\link{get_TRR}} - Identify entry regions
#' }
#'
#' \strong{Downstream Analysis:}
#' \itemize{
#'   \item \code{\link{gene_comparison_tmaker_adj_all}} - Differential expression
#'   \item \code{\link{gene_comparison_umi_adj}} - UMI-adjusted comparison
#' }
#'
#' @section Workflow:
#'
#' A typical analysis follows these steps:
#' \enumerate{
#'   \item Load and prepare data with \code{data_preparation_TRR()}
#'   \item Identify tumor clusters with \code{plot_seurat_clusters()}
#'   \item Estimate J1 with \code{estimate_J1_dat()} and \code{find_optimal_mu_Q()}
#'   \item Decompose spatial data with \code{get_lines()}
#'   \item Detect change points with \code{get_cps_and_save()}
#'   \item Calculate slopes with \code{calculate_slopes()}
#'   \item Estimate J2 with \code{plot_J2()}
#'   \item Run Gibbs sampling with \code{gibbs_sampling_new()}
#'   \item Calculate probabilities with \code{calculate_final_prob()}
#'   \item Identify TRR with \code{get_TRR()}
#'   \item Perform downstream analysis
#' }
#'
#' @docType package
#' @name TRRident-package
#' @aliases TRRident
#'
#' @useDynLib TRRident, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
