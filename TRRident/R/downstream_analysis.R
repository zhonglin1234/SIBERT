#' @title Downstream Analysis Functions
#' @description Functions for differential expression analysis between
#' T cell entry regions and other tumor areas
#' @name downstream_analysis
NULL

#' Gene Comparison with T Cell Marker Adjustment
#'
#' Performs differential expression analysis between TRR spots and control spots,
#' adjusting for T cell marker expression and optionally UMI counts using
#' negative binomial regression.
#'
#' @param genes Character vector of gene names to test
#' @param ids0 Character vector of control spot IDs
#' @param ids1 Character vector of TRR spot IDs
#' @param umi.adj Logical; whether to adjust for UMI counts (default: TRUE)
#'
#' @return Data frame with columns:
#' \itemize{
#'   \item gene - Gene name
#'   \item logFC - Log2 fold change
#'   \item p_value - Raw p-value
#'   \item adj_p_value - BH-adjusted p-value
#'   \item log_adj_p_value - -log10 of adjusted p-value
#' }
#'
#' @importFrom MASS glm.nb
#' @export
#'
#' @examples
#' \dontrun{
#' results <- gene_comparison_tmaker_adj_all(
#'   genes = vgg,
#'   ids0 = control_spots,
#'   ids1 = trr_spots,
#'   umi.adj = TRUE
#' )
#' }
gene_comparison_tmaker_adj_all <- function(genes, ids0, ids1, umi.adj = TRUE) {
  # Get required variables from parent scope
  sp.counts.r <- get("sp.counts.r", envir = parent.frame())
  meta_dat <- get("meta_dat", envir = parent.frame())
  tcell_markers <- get("tcell_markers", envir = parent.frame())

  results <- data.frame(
    gene = genes,
    logFC = numeric(length(genes)),
    p_value = numeric(length(genes))
  )
  rownames(results) <- genes

  xx <- cbind(
    c(rep(0, length(ids0)), rep(1, length(ids1))),
    meta_dat[c(ids0, ids1), 'nCount_Spatial'],
    rowSums(sp.counts.r[c(ids0, ids1), tcell_markers])
  )

  for (i in genes) {
    YY <- sp.counts.r[c(ids0, ids1), i]
    tmp.dat <- data.frame(trr = xx[, 1], umi = xx[, 2], tcell_sum = xx[, 3], y = YY)

    tryCatch({
      if (umi.adj) {
        mm <- MASS::glm.nb(y ~ tcell_sum + umi + trr, data = tmp.dat)
      } else {
        mm <- MASS::glm.nb(y ~ tcell_sum + trr, data = tmp.dat)
      }

      beta <- summary(mm)$coef['trr', 'Estimate']
      logFC <- log2(exp(beta))
      p_value <- summary(mm)$coef['trr', 'Pr(>|z|)']

      results[i, "logFC"] <- logFC
      results[i, "p_value"] <- p_value
    }, error = function(e) {
      results[i, "logFC"] <- NA
      results[i, "p_value"] <- NA
    })
  }

  results$adj_p_value <- stats::p.adjust(results$p_value, method = "BH")
  results$log_adj_p_value <- -log10(results$adj_p_value)
  results$adj_p_value <- round(results$adj_p_value, 4)

  return(results)
}


#' Gene Comparison with UMI Adjustment Only
#'
#' Performs differential expression analysis using t-tests on normalized
#' count data.
#'
#' @param genes Character vector of gene names to test
#' @param ids0 Character vector of control spot IDs
#' @param ids1 Character vector of TRR spot IDs
#'
#' @return Data frame with differential expression results
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- gene_comparison_umi_adj(
#'   genes = vgg,
#'   ids0 = control_spots,
#'   ids1 = trr_spots
#' )
#' }
gene_comparison_umi_adj <- function(genes, ids0, ids1) {
  sp.counts.norm <- get("sp.counts.norm", envir = parent.frame())

  group_0 <- sp.counts.norm[ids0, genes]
  group_1 <- sp.counts.norm[ids1, genes]

  results <- data.frame(
    gene = genes,
    logFC = numeric(length(genes)),
    p_value = numeric(length(genes))
  )

  for (j in seq_along(genes)) {
    gene_0 <- group_0[, j]
    gene_1 <- group_1[, j]
    message(j)

    ttest <- stats::t.test(gene_0, gene_1)
    m0 <- ttest$estimate[1]
    m1 <- ttest$estimate[2]
    logFC <- log2(m1 / m0)

    results$logFC[j] <- logFC
    results$p_value[j] <- ttest$p.value
  }

  results$adj_p_value <- stats::p.adjust(results$p_value, method = "BH")
  results$log_adj_p_value <- -log10(results$adj_p_value)
  results$adj_p_value_round <- round(results$adj_p_value, 4)

  return(results)
}


#' Compare TRR with Tumor Region
#'
#' Compares gene expression between TRR spots and other tumor spots
#' using t-tests on normalized data.
#'
#' @param sp.dat Expression matrix (normalized)
#' @param tumor.loc Data frame with tumor spot coordinates
#' @param cut.off Log2 fold change cutoff for filtering (default: -1)
#'
#' @return Data frame of significantly downregulated genes in TRR
#'
#' @export
compare_with_tumor <- function(sp.dat, tumor.loc, cut.off = -1) {
  ids1 <- get("ids1", envir = parent.frame())

  dat <- sp.dat[rownames(tumor.loc), ]
  TRR <- rep(0, nrow(dat))
  dat <- cbind(dat, TRR)
  dat[which(rownames(tumor.loc) %in% ids1), 'TRR'] <- 1

  tmp.mat <- matrix(NA, nrow = ncol(sp.dat), ncol = 3)

  for (i in seq_len(ncol(sp.dat))) {
    ttest <- stats::t.test(
      dat[which(rownames(tumor.loc) %in% ids1), colnames(sp.dat)[i]],
      dat[which(!rownames(tumor.loc) %in% ids1), colnames(sp.dat)[i]]
    )
    m1 <- ttest$estimate[1]
    m0 <- ttest$estimate[2]
    tmp.mat[i, ] <- c(m1, m0, ttest$p.value)
  }

  colnames(tmp.mat) <- c('mean_TRR', 'mean_other', "pvalue")
  rownames(tmp.mat) <- colnames(sp.dat)
  tmp.mat <- tmp.mat[order(tmp.mat[, 'pvalue']), ]
  tmp.mat <- data.frame(tmp.mat)
  tmp.mat$adj_p_value <- stats::p.adjust(tmp.mat$pvalue, method = 'BH')
  tmp.mat$logFC <- log2(tmp.mat$mean_TRR / tmp.mat$mean_other)

  low.mat <- tmp.mat[which(tmp.mat$logFC < cut.off & tmp.mat$adj_p_value < 0.05), ]

  return(low.mat)
}
