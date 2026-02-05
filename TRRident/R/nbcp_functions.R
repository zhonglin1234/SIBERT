#' @title Change Point Detection with Negative Binomial Regression
#' @description Functions for detecting change points using PELT algorithm
#' with negative binomial regression and Polya-Gamma augmentation
#' @name nbcp_functions
NULL

#' Optimal Partitioning with PELT
#'
#' Detects change points using the Pruned Exact Linear Time (PELT) algorithm
#' with negative binomial regression.
#'
#' @param n_iter Number of MCMC iterations for parameter estimation
#' @param Y Response variable (count data)
#' @param X Design matrix
#' @param method Method for estimation ('PG' for Polya-Gamma)
#' @param alpha Penalty parameter
#' @param K Additional pruning constant (default: 0)
#'
#' @return A list containing:
#' \itemize{
#'   \item total_cost - Total cost of optimal segmentation
#'   \item segmentation - List of segment start/end indices
#'   \item change_points - Vector of change point locations
#'   \item cost_mat - Cost matrix
#' }
#'
#' @keywords internal
optimal_partitioning_pelt <- function(n_iter = 20, Y, X, method, alpha = 1, K = 0) {
  n <- length(Y)

  cost <- rep(Inf, n + 1)
  partition <- rep(0, n + 1)
  cost[1] <- 0

  R <- list()
  R[[1]] <- 0
  R[[2]] <- 1
  R[[3]] <- 1

  R[[4]] <- c(1, 3)
  R[[5]] <- c(1, 3, 4)

  cost.mat <- matrix(NA, n, n)

  for (end in 2:n) {
    for (start in R[[end]]) {
      if (start == 2) next
      if (end - start + 1 >= 2) {
        current_cost <- segment_cost(n_iter = n_iter, Y = Y, X = X,
                                      start = start, end = end,
                                      method = method, alpha = alpha)
        current_cost_aic <- current_cost[1]
        current_cost_llh <- current_cost[2]
        cost.mat[start, end] <- current_cost_llh

        if (cost[start] + current_cost_aic < cost[end + 1]) {
          cost[end + 1] <- cost[start] + current_cost_aic
          partition[end + 1] <- start
        }
      }
    }

    # Pruning
    new_R <- c()
    if (end >= 5 & end < n) {
      tmp.cp <- R[[end]]
      for (start in tmp.cp[2:length(tmp.cp)]) {
        pruning_cost <- cost.mat[start, end]
        f_t <- cost[start]

        if (f_t + pruning_cost + K < cost[end + 1]) {
          new_R <- c(new_R, start)
        }
      }
      R[[end + 1]] <- unique(c(1, new_R, end))
    }
  }

  # Reconstruct segmentation
  segmentation <- list()
  change_points <- c()
  end <- n
  while (end > 0) {
    start <- partition[end + 1]
    segmentation <- append(segmentation, list(c(start, end)))
    if (start != 1) {
      change_points <- c(change_points, start)
    }
    end <- start - 1
  }

  return(list(
    total_cost = cost[n + 1],
    segmentation = rev(segmentation),
    change_points = rev(change_points),
    cost_mat = cost.mat
  ))
}


#' Segment Cost Function
#'
#' Calculates the cost of a segment using negative binomial regression.
#'
#' @param n_iter Number of MCMC iterations
#' @param Y Response variable
#' @param X Design matrix
#' @param start Start index of segment
#' @param end End index of segment
#' @param method Estimation method
#' @param alpha Penalty parameter
#'
#' @return Vector of (AIC-like cost, log-likelihood cost)
#'
#' @keywords internal
segment_cost <- function(n_iter, Y, X, start, end, method, alpha) {
  Y_seg <- Y[start:end]
  X_seg <- X[start:end, , drop = FALSE]

  # Update X_seg to have sequential position
  X_seg[, 2] <- 1:nrow(X_seg)

  np <- ncol(X_seg) + 1  # Number of parameters (beta + dispersion)

  result <- get_NBR_ests(n_iter = n_iter, Y = Y_seg, X = X_seg, method = method)

  # Calculate negative log-likelihood
  neg_llh <- result$neg_llh

  # AIC-like cost with penalty
  cost_aic <- neg_llh + np * 2 * alpha

  return(c(cost_aic, neg_llh))
}


#' Get Negative Binomial Regression Estimates
#'
#' Estimates parameters of negative binomial regression using
#' Polya-Gamma augmentation for Bayesian inference.
#'
#' @param n_iter Number of MCMC iterations
#' @param Y Response variable (counts)
#' @param X Design matrix
#' @param method Estimation method ('PG' for Polya-Gamma)
#'
#' @return List with beta coefficients and negative log-likelihood
#'
#' @keywords internal
get_NBR_ests <- function(n_iter, Y, X, method = 'PG') {
  n <- length(Y)
  p <- ncol(X)

  # Prior parameters
  b <- rep(0, p)
  B <- diag(100, p)

  # Initial values
  beta <- rep(0, p)
  phi <- 1

  # Run Polya-Gamma sampler
  if (method == 'PG') {
    beta_post <- cpp_update_beta(
      XX = X,
      y_obs = Y,
      b = b,
      B = B,
      beta = beta,
      phi = phi,
      niter_beta = n_iter
    )

    phi_post <- cpp_update_phi_MCMC(
      XX = X,
      y_obs = Y,
      beta = beta_post,
      a_phi = 0.01,
      b_phi = 0.01,
      sd_phi = 0.5,
      phi = phi,
      niter_phi = n_iter
    )

    # Calculate negative log-likelihood at posterior means
    psi <- X %*% beta_post
    pp <- 1 / (1 + exp(-psi))
    lambda <- phi_post * pp / (1 - pp)

    neg_llh <- -sum(lgamma(Y + phi_post) - lgamma(Y + 1) - lgamma(phi_post) +
                      phi_post * (log(phi_post) - log(lambda + phi_post)) +
                      Y * (log(lambda) - log(lambda + phi_post)))

    return(list(beta = beta_post, phi = phi_post, neg_llh = neg_llh))
  }

  stop("Unknown method: ", method)
}
