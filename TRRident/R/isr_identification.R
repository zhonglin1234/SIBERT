#' @title T Cell Entry Region Identification
#' @description Functions for identifying T cell entry regions (TRR) using
#' Gibbs sampling on an Ising model with impurities
#' @name isr_identification
NULL

#' Gibbs Sampling for TRR Identification
#'
#' Performs Gibbs sampling on the Ising model with impurities to estimate
#' posterior probabilities of each spot belonging to a T cell entry region.
#'
#' @param ss Initial spin configuration (vector of -1 and 1)
#' @param spot_names Character vector of spot names
#' @param n_add Named numeric vector of additional connections at impurity sites
#' @param nb_list Named list of neighbors for each spot
#' @param new_nb_list Modified neighbor list for impurity spots
#' @param J Coupling parameter (J1)
#' @param loc Data frame with spot coordinates
#' @param iterations Number of Gibbs sampling iterations
#'
#' @return Matrix of spin configurations (iterations x spots)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ss.mat <- gibbs_sampling_new(
#'   ss = ss,
#'   spot_names = names(obs_ss),
#'   n_add = n.add * J2,
#'   nb_list = nb.r,
#'   new_nb_list = new.nb.list,
#'   J = J1,
#'   loc = loc.raw,
#'   iterations = 5000
#' )
#' }
gibbs_sampling_new <- function(ss, spot_names, n_add, nb_list, new_nb_list,
                                J, loc, iterations) {
  names(ss) <- spot_names
  n <- length(ss)
  ss_samples <- matrix(0, nrow = iterations, ncol = n)
  colnames(ss_samples) <- spot_names
  ss_samples[1, ] <- ss

  for (iter in 2:iterations) {
    if (iter %% 100 == 0) message(iter)
    current_ss <- ss_samples[iter - 1, ]

    for (i in 1:n) {
      prob_1 <- get_cond_pp_new(
        spot.name = spot_names[i],
        n.add = n_add,
        nb.list = nb_list,
        new.nb.list = new_nb_list,
        ss = current_ss,
        J = J
      )
      current_ss[i] <- sample(c(1, -1), size = 1, prob = c(prob_1, 1 - prob_1))
    }
    ss_samples[iter, ] <- current_ss
  }

  return(ss_samples)
}


#' Calculate Final Probabilities
#'
#' Calculates posterior probabilities of belonging to TRR for each spot
#' from Gibbs sampling output, after discarding burn-in samples.
#'
#' @param ss_samples Matrix of spin configurations from gibbs_sampling_new()
#' @param loc Data frame with spot coordinates
#' @param burn_in Number of burn-in iterations to discard
#' @param J Coupling parameter (for plot title)
#'
#' @return Named numeric vector of posterior probabilities
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ss.mat[ss.mat == -1] <- 0
#' final.probs <- calculate_final_prob(ss.mat, loc = loc.dat, burn_in = 2000, J = J1)
#' }
calculate_final_prob <- function(ss_samples, loc, burn_in, J) {
  post_burn_in_samples <- ss_samples[(burn_in + 1):nrow(ss_samples), ]
  final_prob <- colMeans(post_burn_in_samples)

  graphics::plot(loc$x, loc$y,
                 col = grDevices::rgb(final_prob, 0, 1 - final_prob),
                 pch = 16, cex = 0.5,
                 xlim = c(min(loc[, 1]) - 50, max(loc[, 1]) + 10000),
                 xlab = "X Coordinate", ylab = "Y Coordinate",
                 main = paste0("J = ", J))

  legend_x <- max(loc$x) + 1000
  legend_y <- max(loc$y)

  graphics::legend(legend_x, legend_y,
                   legend = c("Low", "High"),
                   fill = c("blue", "red"),
                   title = "Pr(ISR=1)", bty = "n", cex = 0.7, x.intersp = 0.5)

  return(final_prob)
}


#' Get T Cell Entry Regions
#'
#' Identifies connected regions of spots with high probability of being
#' T cell entry regions, filtering by minimum size and tumor location.
#'
#' @param final.probs Named numeric vector of posterior probabilities
#'
#' @return List of character vectors, each containing spot names for a TRR
#'
#' @importFrom igraph graph_from_adjacency_matrix components
#' @export
#'
#' @examples
#' \dontrun{
#' final.regions <- get_TRR(final.probs)
#' }
get_TRR <- function(final.probs) {
  # Get required variables from parent scope
  loc.raw <- get("loc.raw", envir = parent.frame())
  r.unit <- get("r.unit", envir = parent.frame())
  tumor.loc <- get("tumor.loc", envir = parent.frame())

  cutp <- stats::quantile(final.probs, probs = seq(0, 1, length = 11))[10]
  ltsa <- final.probs
  ltsa[final.probs > cutp] <- 1
  ltsa[final.probs <= cutp] <- 0

  red_spots <- loc.raw[which(rownames(loc.raw) %in% names(ltsa)[which(ltsa == 1)]), ]

  # Create graph for connected components
  threshold <- 1.5 * r.unit
  dist_matrix <- as.matrix(stats::dist(red_spots))
  adj_matrix <- dist_matrix <= threshold & dist_matrix > 0
  graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

  # Identify connected components
  components <- igraph::components(graph)

  # Filter regions with >= 10 spots
  regions <- lapply(1:length(components$csize), function(x)
    names(components$membership)[which(components$membership == x)])
  final.regions <- regions[which(components$csize >= 10)]

  # Filter to tumor region
  final.regions <- lapply(final.regions, function(x)
    x[which(x %in% rownames(tumor.loc))])
  tmp.ll <- unlist(lapply(final.regions, length))
  final.regions <- final.regions[which(tmp.ll >= 10)]

  # Visualization
  graphics::par(mfrow = c(1, 2))

  graphics::plot(loc.raw, main = "Spatial Plot", xlab = "X", ylab = "Y",
                 pch = 16, col = "gray", ylim = c(0, max(loc.raw[, 2] + 2000)))
  graphics::points(tumor.loc, col = 'yellow', cex = 0.6)
  graphics::points(loc.raw[which(rownames(loc.raw) %in% names(ltsa)[which(ltsa == 1)]), ],
                   col = "red", pch = 19, cex = 0.5)
  graphics::legend("topright", legend = c("LTSA", "Tumor Area"),
                   col = c("red", "yellow"), pch = 19, pt.cex = 1.2, bty = "n")

  graphics::plot(loc.raw, main = "Filtered TRR Regions", xlab = "X", ylab = "Y",
                 col = "gray", pch = 16, ylim = c(0, max(loc.raw[, 2] + 2000)))
  graphics::points(tumor.loc, col = 'yellow', cex = 0.6)
  graphics::points(loc.raw[unlist(final.regions), ], col = "red", pch = 19, cex = 0.7)
  graphics::legend("topright", legend = c("Entry region", "Tumor Area"),
                   col = c("red", "yellow"), pch = 19, pt.cex = 1.2, bty = "n")

  graphics::par(mfrow = c(1, 1))

  return(final.regions)
}


# Internal helper functions

#' Get Endpoints from Segments
#' @keywords internal
get_endpoints <- function(all_segs, spots) {
  endpoints <- matrix(-1, nrow = length(spots), ncol = length(all_segs))
  rownames(endpoints) <- spots
  colnames(endpoints) <- c('degree_30', 'degree_90', 'degree_150')

  for (i in seq_along(all_segs)) {
    tmp.segs <- all_segs[[i]]
    tmp.endpoints <- unlist(lapply(tmp.segs, function(x)
      unlist(lapply(x, function(y) rownames(y)[c(1, nrow(y))]))))
    endpoints[tmp.endpoints, i] <- 1
  }

  endpoints <- cbind(endpoints, endpoints)
  colnames(endpoints)[4:6] <- c('degree_-150', 'degree_-90', 'degree_-30')

  return(endpoints)
}


#' Calculate Direction from Center to Neighbor
#' @keywords internal
calculate_direction_center_to_nb <- function(spot_i, neighbors, loc) {
  spot_i_coords <- loc[spot_i, ]
  neighbor_coords <- loc[neighbors, ]

  dx <- neighbor_coords[, 1] - as.numeric(spot_i_coords[1])
  dy <- neighbor_coords[, 2] - as.numeric(spot_i_coords[2])

  angle <- atan2(dy, dx) * 180 / pi
  angle <- round(angle, -1)
  angle <- paste0('degree_', angle)

  return(angle)
}


#' Get Endpoints Neighbors
#' @keywords internal
get_endpoints_neighbors <- function(endpoints.mat, slope.dat, nb.list, loc, margin.spots) {
  nb.of.eps.list <- list()
  endp.works.list <- list()
  angles <- colnames(endpoints.mat)

  for (i in 1:6) {
    tmp.eps <- rownames(endpoints.mat)[which(endpoints.mat[, i] == 1)]
    tmp.eps <- tmp.eps[which(!tmp.eps %in% margin.spots)]
    aa <- angles[i]
    nb_of_ends <- c()
    endp_works <- c()

    for (j in tmp.eps) {
      nb <- nb.list[[j]]
      nb <- nb[which(!nb %in% tmp.eps)]
      direction <- calculate_direction_center_to_nb(j, nb, loc)

      if (aa %in% direction) {
        if (slope.dat[nb[which(direction == aa)], aa] == 1) {
          endp_works <- c(endp_works, j)
          nb_of_ends <- c(nb_of_ends, nb[which(direction == aa)])
        }
      }
    }

    names(nb_of_ends) <- endp_works
    endp.works.list <- c(endp.works.list, list(endp_works))
    nb.of.eps.list <- c(nb.of.eps.list, list(nb_of_ends))
  }

  names(nb.of.eps.list) <- colnames(endpoints.mat)
  names(endp.works.list) <- colnames(endpoints.mat)

  all_spots <- unlist(nb.of.eps.list, use.names = TRUE)
  adj_ep_spots <- split(names(all_spots), all_spots)
  adj_ep_spots_noangle <- lapply(adj_ep_spots, function(x)
    substr(x, nchar(x) - 17, nchar(x)))
  spots.add <- unlist(lapply(adj_ep_spots_noangle, length))

  new.nb.list <- list()
  for (spot.name in names(spots.add)) {
    nb <- nb.list[[spot.name]]
    new.nb <- nb[which(!nb %in% adj_ep_spots_noangle[[spot.name]])]
    new.nb.list <- c(new.nb.list, list(new.nb))
  }

  return(list(n.add = spots.add, new.nb.list = new.nb.list))
}


#' Get Conditional Probability for Gibbs Sampling
#' @keywords internal
get_cond_pp_new <- function(spot.name, n.add, nb.list, new.nb.list, ss, J) {
  if (spot.name %in% names(n.add)) {
    nb <- new.nb.list[[spot.name]]
    ss.nb <- ss[nb]
    Ham1 <- -J * (length(which(ss.nb == 1)) + n.add[spot.name])
    Ham0 <- -J * (length(which(ss.nb == -1)))
    pp1 <- exp(-Ham1) / (exp(-Ham1) + exp(-Ham0))
  } else {
    nb <- nb.list[[spot.name]]
    ss.nb <- ss[nb]
    Ham1 <- -J * (length(which(ss.nb == 1)))
    Ham0 <- -J * (length(which(ss.nb == -1)))
    pp1 <- exp(-Ham1) / (exp(-Ham1) + exp(-Ham0))
  }

  return(pp1)
}
