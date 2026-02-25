#' @title Line Management and Change Point Detection
#' @description Functions for decomposing 2D spatial data into 1D lines
#' and detecting change points in gene expression
#' @name line_management
NULL

#' Get Lines from Spatial Data
#'
#' Decomposes 2D spatial transcriptomics data into 1D lines along multiple
#' directions. For hexagonal grids (Visium), extracts lines along 30°, 90°,
#' and 150° directions.
#'
#' @param loc Data frame with x and y coordinates
#' @param grid.type Type of spatial grid: 'square' or 'hex'
#' @param aa Tolerance parameter for grouping coordinates (default: 0.2)
#' @param max_gap Maximum allowed gap between adjacent spots on a line,
#'   as a multiple of minimum spot distance (default: 10)
#'
#' @return A list of line lists:
#' \itemize{
#'   \item For 'hex': lines_30, lines_90, lines_150
#'   \item For 'square': vertical_lines, horizontal_lines
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' all.lines <- get_lines(loc.raw, grid.type = 'hex', max_gap = 7, aa = 0.2)
#' plot_lines(loc.raw, all.lines, grid.type = 'hex')
#' }
get_lines <- function(loc, grid.type, aa = 0.2, max_gap = 10) {
  loc <- data.frame(loc[, 1:2])
  dist.mat <- as.matrix(stats::dist(loc))
  min.dist <- min(dist.mat[upper.tri(dist.mat)])

  # Function to split lines by max_gap
  split_by_gap <- function(line_data, dist_cutoff) {
    if (nrow(line_data) < 2) return(list())

    segments <- list()
    current_segment <- line_data[1, , drop = FALSE]

    for (i in 2:nrow(line_data)) {
      point_dist <- sqrt(sum((line_data[i, ] - line_data[i - 1, ])^2))
      if (point_dist > dist_cutoff) {
        if (nrow(current_segment) >= 2) {
          segments <- append(segments, list(current_segment))
        }
        current_segment <- line_data[i, , drop = FALSE]
      } else {
        current_segment <- rbind(current_segment, line_data[i, , drop = FALSE])
      }
    }

    if (nrow(current_segment) >= 2) {
      segments <- append(segments, list(current_segment))
    }

    return(segments)
  }

  if (grid.type == 'square') {
    loc$vertical_group <- group_close_coords_square(loc[, 1], cutoff = min.dist * aa)
    loc$horizontal_group <- group_close_coords_square(loc[, 2], cutoff = min.dist * aa)

    vertical_lines <- split(loc[, 1:2], loc$vertical_group)
    vertical_lines <- lapply(vertical_lines, function(df) df[order(df[, 2]), ])
    vertical_lines <- unlist(lapply(vertical_lines, function(line)
      split_by_gap(line, max_gap * min.dist)), recursive = FALSE)
    vertical_lines <- vertical_lines[sapply(vertical_lines, nrow) >= 2]
    names(vertical_lines) <- paste0("Vertical Line ", seq_along(vertical_lines))

    horizontal_lines <- split(loc[, 1:2], loc$horizontal_group)
    horizontal_lines <- lapply(horizontal_lines, function(df) df[order(df[, 1]), ])
    horizontal_lines <- unlist(lapply(horizontal_lines, function(line)
      split_by_gap(line, max_gap * min.dist)), recursive = FALSE)
    horizontal_lines <- horizontal_lines[sapply(horizontal_lines, nrow) >= 2]
    names(horizontal_lines) <- paste0("Horizontal Line ", seq_along(horizontal_lines))

    out.lines <- list(vertical_lines = vertical_lines, horizontal_lines = horizontal_lines)
    return(out.lines)
  }

  if (grid.type == 'hex') {
    rotate_30 <- matrix(c(cos(pi / 3), -sin(pi / 3), sin(pi / 3), cos(pi / 3)), ncol = 2)
    rotate_150 <- matrix(c(cos(2 * pi / 3), -sin(2 * pi / 3),
                           sin(2 * pi / 3), cos(2 * pi / 3)), ncol = 2)

    loc_rotated_30 <- as.matrix(loc) %*% rotate_30
    loc_rotated_150 <- as.matrix(loc) %*% rotate_150

    loc$group_30 <- group_close_coords_hex(loc_rotated_30[, 1], cutoff = min.dist * aa)
    loc$group_150 <- group_close_coords_hex(loc_rotated_150[, 1], cutoff = min.dist * aa)
    loc$group_90 <- group_close_coords_hex(loc[, 1], cutoff = min.dist * aa)

    lines_30 <- split(loc[, 1:2], loc$group_30)
    lines_30 <- lapply(lines_30, function(df) df[order(df[, 2]), ])
    lines_30 <- unlist(lapply(lines_30, function(line)
      split_by_gap(line, max_gap * min.dist)), recursive = FALSE)
    lines_30 <- lines_30[sapply(lines_30, nrow) >= 2]

    lines_150 <- split(loc[, 1:2], loc$group_150)
    lines_150 <- lapply(lines_150, function(df) df[order(df[, 2]), ])
    lines_150 <- unlist(lapply(lines_150, function(line)
      split_by_gap(line, max_gap * min.dist)), recursive = FALSE)
    lines_150 <- lines_150[sapply(lines_150, nrow) >= 2]

    lines_90 <- split(loc[, 1:2], loc$group_90)
    lines_90 <- lapply(lines_90, function(df) df[order(df[, 2]), ])
    lines_90 <- unlist(lapply(lines_90, function(line)
      split_by_gap(line, max_gap * min.dist)), recursive = FALSE)
    lines_90 <- lines_90[sapply(lines_90, nrow) >= 2]

    out.lines <- list(lines_30 = lines_30, lines_90 = lines_90, lines_150 = lines_150)
    return(out.lines)
  }
}


#' Plot Lines on Spatial Data
#'
#' Visualizes the decomposed lines overlaid on the spatial spot coordinates.
#'
#' @param loc Data frame with x and y coordinates
#' @param lines.list Output from get_lines()
#' @param grid.type Type of spatial grid: 'square' or 'hex'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' all.lines <- get_lines(loc.raw, grid.type = 'hex')
#' plot_lines(loc.raw, all.lines, grid.type = 'hex')
#' }
plot_lines <- function(loc, lines.list, grid.type) {
  if (grid.type == 'square') {
    graphics::par(mfrow = c(1, 2))

    vertical_lines <- lines.list$vertical_lines
    graphics::plot(loc[, 1:2], xlab = "X Coordinate", ylab = "Y Coordinate",
                   pch = 19, cex = 0.3, col = 'grey', main = 'Vertical lines')
    for (i in seq_along(vertical_lines)) {
      graphics::lines(vertical_lines[[i]][, 1], vertical_lines[[i]][, 2],
                      col = "blue", lwd = 0.5)
    }

    horizontal_lines <- lines.list$horizontal_lines
    graphics::plot(loc[, 1:2], xlab = "X Coordinate", ylab = "Y Coordinate",
                   pch = 19, cex = 0.3, col = 'grey', main = 'Horizontal lines')
    for (i in seq_along(horizontal_lines)) {
      graphics::lines(horizontal_lines[[i]][, 1], horizontal_lines[[i]][, 2],
                      col = "red", lwd = 0.5)
    }
  }

  if (grid.type == 'hex') {
    graphics::par(mfrow = c(1, 3))

    lines_30 <- lines.list$lines_30
    graphics::plot(loc[, 1:2], xlab = "X Coordinate", ylab = "Y Coordinate",
                   pch = 19, cex = 0.3, col = 'grey', main = '30 degree lines')
    for (i in seq_along(lines_30)) {
      graphics::lines(lines_30[[i]][, 1], lines_30[[i]][, 2], col = "blue", lwd = 0.5)
    }

    lines_90 <- lines.list$lines_90
    graphics::plot(loc[, 1:2], xlab = "X Coordinate", ylab = "Y Coordinate",
                   pch = 19, cex = 0.3, col = 'grey', main = '90 degree lines')
    for (i in seq_along(lines_90)) {
      graphics::lines(lines_90[[i]][, 1], lines_90[[i]][, 2], col = "green", lwd = 0.5)
    }

    lines_150 <- lines.list$lines_150
    graphics::plot(loc[, 1:2], xlab = "X Coordinate", ylab = "Y Coordinate",
                   pch = 19, cex = 0.3, col = 'grey', main = '150 degree lines')
    for (i in seq_along(lines_150)) {
      graphics::lines(lines_150[[i]][, 1], lines_150[[i]][, 2], col = "red", lwd = 0.5)
    }
  }
}


#' Get Change Points and Save Results
#'
#' Detects change points in gene expression along each line using
#' PELT algorithm with negative binomial regression.
#'
#' @param all.lines Output from get_lines()
#' @param YY Named numeric vector of response values (e.g., exhaustion marker sum)
#' @param XX Design matrix for regression
#' @param K Pruning parameter for PELT (default: 0)
#' @param alpha_vec Vector of penalty values to test
#' @param method_cp Method for change point detection ('PG' for Polya-Gamma)
#' @param dat_name Base name for output files
#' @param grid.type Type of spatial grid: 'square' or 'hex'
#'
#' @return A list containing change points and segments for each penalty value
#'
#' @export
#'
#' @examples
#' \dontrun{
#' all.output <- get_cps_and_save(
#'   all.lines = all.lines,
#'   YY = YY,
#'   XX = XX,
#'   K = 0,
#'   alpha_vec = seq(0.1, 1, length = 10),
#'   method_cp = 'PG',
#'   dat_name = 'data/melanoma1_tcelladj',
#'   grid.type = 'hex'
#' )
#' }
get_cps_and_save <- function(all.lines, YY, XX, K = 0, alpha_vec,
                              method_cp = 'PG', dat_name, grid.type) {
  all.cps.list <- list()
  all.segs.list <- list()

  for (a in seq_along(alpha_vec)) {
    alpha <- alpha_vec[a]
    message(paste("Processing penalty alpha =", alpha))

    all.cps <- list()
    all.segs <- list()

    for (i in seq_along(all.lines)) {
      tmp.lines <- all.lines[[i]]
      tmp.cps <- list()
      tmp.segs <- list()

      for (j in seq_along(tmp.lines)) {
        line_spots <- rownames(tmp.lines[[j]])

        if (length(line_spots) < 4) {
          tmp.cps[[j]] <- integer(0)
          tmp.segs[[j]] <- list(tmp.lines[[j]])
          next
        }

        # Prepare data for change point detection
        Y_line <- YY[line_spots]
        X_line <- XX[line_spots, ]
        X_line[, 2] <- 1:length(line_spots)

        # Run PELT with negative binomial regression
        result <- tryCatch({
          optimal_partitioning_pelt(
            n_iter = 20,
            Y = Y_line,
            X = X_line,
            method = method_cp,
            alpha = alpha,
            K = K
          )
        }, error = function(e) {
          list(change_points = integer(0), segmentation = list(c(1, length(line_spots))))
        })

        tmp.cps[[j]] <- result$change_points
        # Convert segmentation to actual data frames
        segs <- lapply(result$segmentation, function(seg) {
          tmp.lines[[j]][seg[1]:seg[2], , drop = FALSE]
        })
        tmp.segs[[j]] <- segs
      }

      all.cps[[i]] <- tmp.cps
      all.segs[[i]] <- tmp.segs
    }

    all.cps.list[[a]] <- all.cps
    all.segs.list[[a]] <- all.segs
  }

  # Save results
  alpha_str <- paste(sprintf("%02d", as.integer(alpha_vec * 100)), collapse = "_")
  save_file <- paste0(dat_name, "_exhst_cps_penalty_", alpha_str, ".Rdata")
  save(all.cps.list, all.segs.list, alpha_vec, file = save_file)

  message(paste("Results saved to:", save_file))

  return(list(all.cps.list = all.cps.list, all.segs.list = all.segs.list))
}


#' Calculate Slopes for Each Segment
#'
#' Calculates the slope of exhaustion marker expression for each segment
#' identified by change point detection.
#'
#' @param all_segs Segments from get_cps_and_save()
#' @param all.lines Lines from get_lines()
#' @param y_markers Character vector of marker genes for response variable
#' @param sp.counts.r Raw count matrix
#' @param method_beta Method for slope estimation ('PG' for Polya-Gamma)
#' @param alpha Penalty value used
#' @param dat_name Base name for output files
#'
#' @return A list with slope.mat and mean.mat matrices
#'
#' @export
calculate_slopes <- function(all_segs, all.lines, y_markers, sp.counts.r,
                              method_beta = 'PG', alpha, dat_name) {
  loc.raw <- get("loc.raw", envir = parent.frame())
  tcell_markers <- get("tcell_markers", envir = parent.frame())
  ex_markers <- get("ex_markers", envir = parent.frame())

  # Get adjusted Y values
  YY <- rowSums(sp.counts.r[, y_markers])
  NN <- length(YY)
  XX <- cbind(rep(1, NN), rep(NA, NN), rowSums(sp.counts.r[, tcell_markers]))
  rownames(XX) <- names(YY)
  dat <- data.frame(x = XX[, 3], y = YY)
  mm <- MASS::glm.nb(y ~ x, data = dat)
  yy_adjusted <- stats::residuals(mm, type = 'pearson')
  names(yy_adjusted) <- names(YY)

  # Calculate slopes for each segment
  all_betas <- list()
  all_means <- list()

  for (i in seq_along(all_segs)) {
    tmp.segs <- all_segs[[i]]

    tmp.betas <- lapply(tmp.segs, function(line_segs) {
      lapply(line_segs, function(seg) {
        if (length(unique(YY[rownames(seg)])) == 1) {
          return(0)
        } else if (nrow(seg) < 3) {
          return(NA)
        } else {
          XX_modified <- XX[rownames(seg), ]
          XX_modified[, 2] <- 1:nrow(seg)
          return(get_NBR_ests(n_iter = 50, Y = YY[rownames(seg)],
                              X = XX_modified, method = method_beta)$beta[2])
        }
      })
    })

    tmp.means <- lapply(tmp.segs, function(line_segs) {
      lapply(line_segs, function(seg) {
        return(mean(yy_adjusted[rownames(seg)]))
      })
    })

    all_betas <- c(all_betas, list(tmp.betas))
    all_means <- c(all_means, list(tmp.means))
  }

  # Assign slopes and means to each spot
  slope.mat <- matrix(NA, nrow = nrow(loc.raw), ncol = length(all.lines))
  mean.mat <- matrix(NA, nrow = nrow(loc.raw), ncol = length(all.lines))
  rownames(slope.mat) <- rownames(loc.raw)
  rownames(mean.mat) <- rownames(loc.raw)

  for (i in seq_along(all_segs)) {
    tmp_segs <- all_segs[[i]]
    tmp_betas <- all_betas[[i]]
    tmp_means <- all_means[[i]]

    for (j in seq_along(tmp_segs)) {
      tmp.segs <- tmp_segs[[j]]
      tmp.beta <- tmp_betas[[j]]
      tmp.mean <- tmp_means[[j]]

      l.segs <- unlist(lapply(tmp.segs, nrow))
      if (min(l.segs) < 3) {
        tmp.beta[which(l.segs < 3)] <- NA
        tmp.mean[which(l.segs < 3)] <- NA
      }

      for (k in seq_along(tmp.segs)) {
        slope.mat[rownames(tmp.segs[[k]]), i] <- tmp.beta[k][[1]]
        mean.mat[rownames(tmp.segs[[k]]), i] <- tmp.mean[k][[1]]
      }
    }
  }

  slope.output <- list(slope.mat = slope.mat, mean.mat = mean.mat)

  alpha_str <- sprintf("%02d", as.integer(alpha * 100))
  save_file <- paste0(dat_name, "_exhst_slopes_penalty", alpha_str, ".Rdata")
  save(slope.output, file = save_file)

  return(slope.output)
}


# Internal helper functions

group_close_coords_square <- function(coords, cutoff) {
  original_order <- order(coords)
  sorted_coords <- sort(coords)
  group <- rep(1, length(sorted_coords))

  for (i in 2:length(sorted_coords)) {
    if (abs(sorted_coords[i] - sorted_coords[i - 1]) >= cutoff) {
      group[i] <- group[i - 1] + 1
    } else {
      group[i] <- group[i - 1]
    }
  }

  group_in_original_order <- rep(NA, length(group))
  group_in_original_order[original_order] <- group

  return(group_in_original_order)
}

group_close_coords_hex <- function(coords, cutoff) {
  original_order <- order(coords)
  sorted_coords <- sort(coords)
  group <- rep(1, length(sorted_coords))

  for (i in 2:length(sorted_coords)) {
    if (abs(sorted_coords[i] - sorted_coords[i - 1]) >= cutoff) {
      group[i] <- group[i - 1] + 1
    } else {
      group[i] <- group[i - 1]
    }
  }

  group_in_original_order <- rep(NA, length(group))
  group_in_original_order[original_order] <- group

  return(group_in_original_order)
}
