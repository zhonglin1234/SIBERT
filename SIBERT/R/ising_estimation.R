#' @title Ising Model Parameter Estimation
#' @description Functions for estimating J1 and J2 parameters of the Ising model
#' @name ising_estimation
NULL

#' Estimate J1 Parameter for Ising Model
#'
#' Estimates the first coupling parameter (J1) for the Ising model using
#' spatial expression data. This parameter controls the interaction strength
#' between neighboring spots.
#'
#' @param genes Character vector of gene names to use for estimation (e.g., 'B2M')
#' @param sp.dd Count matrix (spots x genes)
#' @param loc.raw
#' @param r.unit
#' @param adj.umi Logical; whether to adjust for UMI counts (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item dat - Data list with location, observed states, and neighbor information
#'   \item sum_obs - Observed sum of interactions
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' j1_data <- estimate_J1_dat(genes = 'B2M', sp.dd = sp.counts.r, adj.umi = FALSE)
#' }
estimate_J1_dat <- function(genes, sp.dd, loc.raw, r.unit, adj.umi = TRUE) {
  graphics::par(mfrow = c(1, 1))
  hla.c1 <- sp.dd[, genes]
  if (length(genes) > 1) {
    hla.c1 <- rowSums(hla.c1)
  }

  if (adj.umi) {
    hla.c1 <- get_adjust_YY(hla.c1)
  }

  plot_spot_by_colors1(hla.c1, loc = loc.raw, leg = 'Level of antigen')

  dat <- get_nb(xx = hla.c1, loc = loc.raw, threshold = r.unit * 1.5)
  obs_ss <- dat[[3]]
  loc.dat <- dat[[2]]
  nb.r <- dat[[4]]
  nb.info <- dat[[5]]

  message("Number of neighbors per spot:")
  print(table(unlist(lapply(nb.r, length))))

  n.pairs <- sum(unlist(lapply(nb.info, length))) / 2
  cat(paste0("Number of pairs: ", n.pairs, "; "))
  sum_obs <- calculate_sum_cpp(obs_ss, nb.info)
  cat(paste0("Number of pairs with same state: ", sum_obs, '; '))
  omega <- sum_obs / n.pairs
  cat(paste0("omega: ", omega))

  return(list(dat, sum_obs))
}


#' Find Optimal mu_Q for Ising Model
#'
#' Uses variational inference to find the optimal mu_Q parameter that
#' minimizes the ELBO (Evidence Lower Bound) for the Ising model.
#'
#' @param mu_Q_candidates Numeric vector of candidate values for mu_Q
#' @param init_x Initial spin configuration (vector of -1 and 1)
#' @param nb.list.cpp Neighbor list formatted for C++ (0-indexed)
#' @param mu Prior mean for J
#' @param sigma2 Prior variance for J
#' @param sigma2_Q Variance for the variational distribution
#' @param K Number of samples for Monte Carlo estimation
#' @param s_obs Observed sum of neighbor interactions
#' @param n_flip Number of spin flips in MCMC
#' @param M Number of samples for gradient estimation
#'
#' @return Named numeric vector of absolute first derivative values for each mu_Q candidate
#'
#' @export
#'
#' @examples
#' \dontrun{
#' first_dr <- find_optimal_mu_Q(
#'   mu_Q_candidates = seq(0.16, 0.44, by = 0.01),
#'   init_x = sample(c(-1, 1), size = 100, prob = c(0.5, 0.5), replace = TRUE),
#'   nb.list.cpp = nb.info,
#'   mu = 0.3, sigma2 = 0.0049, sigma2_Q = 0.0001,
#'   K = 10, s_obs = sum_obs, n_flip = 60000, M = 10000
#' )
#' J1 <- as.numeric(names(first_dr)[which.min(first_dr)])
#' }
find_optimal_mu_Q <- function(mu_Q_candidates = seq(0.16, 0.44, by = 0.01),
                               init_x, nb.list.cpp, mu, sigma2, sigma2_Q,
                               K, s_obs, n_flip, M) {
  first_derivative_values <- c()

  for (i in seq_along(mu_Q_candidates)) {
    mu_Q <- mu_Q_candidates[i]

    # Sample epsilon_k ~ N(0, 1) and compute J_k
    J_k <- stats::rnorm(K, mean = mu_Q, sd = sqrt(sigma2_Q))

    d_logZ_dJk <- sapply(J_k, function(j) {
      MC_sampling_x_cpp(init_x = init_x, nb_list = nb.list.cpp,
                        J_k = j, n_flip = n_flip, M = M)[[2]]
    })

    # Calculate the first derivative of the ELBO with respect to mu_Q
    first_derivative <- s_obs - (1 / K) * sum(d_logZ_dJk) - (mu_Q - mu) / sigma2

    first_derivative_values <- c(first_derivative_values, abs(first_derivative))
    nn <- length(first_derivative_values)

    # Early stopping if gradient is increasing
    if (nn > 5) {
      if (min(diff(first_derivative_values[(nn - 3):nn])) > 0) {
        break
      }
    }
    print(first_derivative_values)
  }

  names(first_derivative_values) <- mu_Q_candidates[1:length(first_derivative_values)]
  return(first_derivative_values)
}


#' Plot J2 Parameter Estimation
#'
#' Estimates the J2 parameter (impurity coupling) by examining the phase
#' transition behavior of the Ising model with impurities.
#'
#' @param J1 The estimated J1 parameter
#' @param lambda.list Vector of candidate J2 values to test
#'
#' @return A list containing:
#' \itemize{
#'   \item mean.list - Percentage of spots with Pr(TRR) > 0.7 for each J2
#'   \item n.add - Number of additional connections at impurity sites
#'   \item new.nb.list - Modified neighbor list accounting for impurities
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' j2_data <- plot_J2(J1 = 0.38, lambda.list = seq(1, 10, length = 5))
#' }
plot_J2 <- function(J1, lambda.list) {
  # Note: This function requires several objects in scope
  # (all.cps.list, all.segs.list, slope.output, etc.)
  # Will need proper parameter passing in full refactor

  all.cps <- all.cps.list[[10]]
  all_segs <- all.segs.list[[10]]
  spots <- rownames(loc.raw)
  endpoints.mat <- get_endpoints(all_segs, spots = spots)

  slope.mat <- slope.output[[1]]
  message("Matrix of slopes:")
  print(utils::head(slope.mat))

  slope.dat <- cbind(slope.mat, -slope.mat)
  colnames(slope.dat) <- paste0('degree_', c(30, 90, 150, -150, -90, -30))
  slope.dat[is.na(slope.dat)] <- 0
  slope.dat[which(slope.dat > -0.01)] <- -1
  slope.dat[which(slope.dat != -1)] <- 1

  message("Number of slopes < -0.01:")
  print(table(slope.dat))

  tmp.dat <- get_endpoints_neighbors(endpoints.mat, slope.dat, nb.r,
                                      loc = loc.raw, margin.spots = margin.spots)
  n.add <- tmp.dat[[1]]
  new.nb.list <- tmp.dat[[2]]

  # Test different J2 values
  ss <- sample(c(-1, 1), size = length(obs_ss), prob = c(0.5, 0.5), replace = TRUE)

  ss.mat.list <- list()
  for (lambda in lambda.list) {
    message(lambda)
    tmp.ss.mat <- gibbs_sampling_new(ss = ss,
                                     spot_names = names(obs_ss),
                                     n_add = n.add * lambda,
                                     nb_list = nb.r,
                                     new_nb_list = new.nb.list,
                                     J = J1,
                                     loc = loc.raw,
                                     iterations = 1000)
    ss.mat.list <- c(ss.mat.list, list(tmp.ss.mat))
  }

  mean.list <- c()
  for (i in seq_along(ss.mat.list)) {
    ss.mat <- ss.mat.list[[i]]
    ss.mat[ss.mat == -1] <- 0
    tmp.final.probs <- calculate_final_prob(ss.mat, loc = loc.dat, burn_in = 500, J = J1)
    mean.list <- c(mean.list, length(which(tmp.final.probs > 0.7)) / length(ss))
  }

  mean.list <- mean.list * 100

  return(list(mean.list = mean.list, n.add = n.add, new.nb.list = new.nb.list))
}


# Internal helper functions

#' Get Neighbors
#' @keywords internal
get_nb <- function(xx, loc, threshold) {
  dist_matrix <- as.matrix(stats::dist(loc))
  nb.list <- lapply(1:nrow(loc), function(i)
    which(dist_matrix[i, ] <= threshold & dist_matrix[i, ] > 0))
  n.nb <- unlist(lapply(nb.list, length))
  spots.keep <- which(n.nb > 0)
  xx.new <- xx[spots.keep]
  loc.new <- loc[spots.keep, ]

  dist_matrix <- as.matrix(stats::dist(loc.new))
  nb.list <- lapply(1:nrow(loc.new), function(i)
    which(dist_matrix[i, ] <= threshold & dist_matrix[i, ] > 0))
  nb.list.cpp <- lapply(nb.list, function(x) as.numeric(x) - 1)

  nb.list <- lapply(1:nrow(loc.new), function(i)
    rownames(loc.new)[which(dist_matrix[i, ] <= threshold & dist_matrix[i, ] > 0)])
  names(nb.list) <- rownames(loc.new)
  names(nb.list.cpp) <- rownames(loc.new)

  gmm <- mclust::Mclust(xx.new, G = 2)

  message("GMM Summary:")
  print(summary(gmm))
  means <- gmm$parameters$mean
  variances <- gmm$parameters$variance$sigmasq
  print(list(means = means, variances = variances))

  clusters <- gmm$classification
  ss <- clusters
  ss[ss == 2] <- -1
  print(table(ss))

  return(list(xx.new, loc.new, ss, nb.list, nb.list.cpp))
}


#' Adjust Y for UMI
#' @keywords internal
get_adjust_YY <- function(YY, XX = NULL) {
  NN <- length(YY)
  if (is.null(XX)) {
    # Default: adjust for total UMI (requires sp.counts.r in scope)
    stop("XX must be provided for UMI adjustment")
  }
  dat <- data.frame(x = XX, y = YY)
  mm <- MASS::glm.nb(y ~ x, data = dat)
  yy_adjusted <- stats::residuals(mm, type = 'pearson')
  names(yy_adjusted) <- names(YY)
  return(yy_adjusted)
}


#' Plot Spots by Color
#'
#' Creates a scatter plot of spatial spots colored by a continuous value.
#'
#' @param loc Data frame with x and y coordinates
#' @param values Numeric vector of values to color by
#' @param leg_title Legend title
#'
#' @return A ggplot object
#'
#' @export
plot_spot_by_colors <- function(loc, values, leg) {
  p1 <- ggplot2::ggplot(loc, ggplot2::aes(x = x, y = y, color = values)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = '',
                  x = "X Coordinate",
                  y = "Y Coordinate",
                  color = leg) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(p1)
}


#' Plot Spots by Color (Base R version)
#'
#' Creates a scatter plot of spatial spots colored by a continuous value
#' using base R graphics.
#'
#' @param loc Data frame with x and y coordinates
#' @param values Numeric vector of values to color by
#' @param leg_title Legend title
#'
#' @export
plot_spot_by_colors1 <- function(loc, values, leg_title) {
  if (!all(c("x", "y") %in% names(loc))) {
    stop("The 'loc' data frame must contain 'x' and 'y' columns.")
  }

  n_colors <- 100
  color_palette <- grDevices::colorRampPalette(c("blue", "red"))(n_colors)

  values_min <- min(values, na.rm = TRUE)
  values_max <- max(values, na.rm = TRUE)

  if (values_max == values_min) {
    normalized_values <- rep(1, length(values))
  } else {
    normalized_values <- (values - values_min) / (values_max - values_min)
  }

  color_indices <- as.numeric(cut(normalized_values, breaks = n_colors, include.lowest = TRUE))
  colors_assigned <- color_palette[color_indices]

  graphics::plot(loc$x, loc$y,
                 col = colors_assigned,
                 pch = 16, cex = 0.5,
                 xlab = "X Coordinate",
                 ylab = "Y Coordinate",
                 main = leg_title)
}
