#' Distance between Cocktail clusters based on direct co-membership phi
#'
#' @description
#' Computes pairwise distances between Cocktail clusters (internal nodes) by
#' computing the \eqn{\phi} coefficient between their binary plot-membership
#' vectors across all plots.
#'
#' Cluster membership is derived from \code{x$Plot.cluster} as:
#' \code{membership = (Plot.cluster > 0)}. This works for both
#' \code{plot_values = "binary"} and \code{plot_values = "rel_cover"}, because
#' in both cases the stored values are positive for plots meeting the Cocktail
#' membership threshold (m-threshold) of a node.
#'
#' For clusters \eqn{A} and \eqn{B}, similarity is:
#' \deqn{\mathrm{sim}(A,B) = \phi(A,B) \in [-1,1],}
#' where \eqn{\phi(A,B)} is computed from the 2×2 contingency table of plot
#' co-membership (plots in \eqn{A} vs. plots in \eqn{B}).
#'
#' The distance is defined as:
#' \deqn{d(A,B) = 1 - \phi(A,B).}
#'
#' Since \eqn{\phi \in [-1,1]}, distances lie in \eqn{[0,2]}.
#'
#' @param x A \code{"cocktail"} object from \code{\link{cocktail_cluster}},
#'   containing at least \code{Plot.cluster} and \code{Cluster.species}.
#'
#' @param clusters Optional cluster identifiers (internal nodes) to be compared.
#'   Can be numeric node IDs (e.g. \code{c(12, 27)}) or character labels
#'   (e.g. \code{c("c_12","c_27")}).
#'   If \code{NULL} or missing, all internal nodes
#'   \code{1:nrow(x$Cluster.species)} are used.
#'
#' @return A \code{\link[stats]{dist}} object with distances \eqn{d = 1-\phi}
#'   between clusters.
#'
#' @details
#' Let \eqn{g_A} and \eqn{g_B} be the binary membership vectors (length = number
#' of plots) for clusters \eqn{A} and \eqn{B}. The \eqn{\phi} coefficient is
#' computed from the 2×2 table:
#' \itemize{
#'   \item \eqn{a}: plots where \eqn{g_A=1} and \eqn{g_B=1},
#'   \item \eqn{b}: plots where \eqn{g_A=1} and \eqn{g_B=0},
#'   \item \eqn{c}: plots where \eqn{g_A=0} and \eqn{g_B=1},
#'   \item \eqn{d}: plots where \eqn{g_A=0} and \eqn{g_B=0}.
#' }
#' \deqn{
#' \phi \;=\; \frac{ad - bc}{\sqrt{(a+c)(b+d)(a+b)(c+d)}}.
#' }
#' Undefined cases (zero denominator) are set to \eqn{\phi=0}.
#'
#' @importFrom stats as.dist
#' @import Matrix
#' @export
cluster_phi_dist <- function(x, clusters = NULL) {

  ## ---- basic checks -------------------------------------------------------
  if (!is.list(x) || !"Cluster.species" %in% names(x)) {
    stop("`x` must be a Cocktail object with a `Cluster.species` component.")
  }
  if (!"Plot.cluster" %in% names(x) || is.null(x$Plot.cluster)) {
    stop("`x$Plot.cluster` is missing; cannot compute cluster distances.")
  }

  CS <- x$Cluster.species
  PC <- x$Plot.cluster

  if (!is.matrix(CS)) stop("`x$Cluster.species` must be a matrix.")

  # allow base matrices but prefer sparse Matrix objects
  if (!inherits(PC, "Matrix")) {
    PC <- Matrix::Matrix(as.matrix(PC), sparse = TRUE)
  }

  n_nodes <- nrow(CS)
  n_plots <- nrow(PC)

  if (n_plots < 1L) stop("`x$Plot.cluster` has zero plots (rows).")
  if (ncol(PC) < n_nodes) {
    stop("`x$Plot.cluster` must have at least nrow(x$Cluster.species) columns (one per node).")
  }

  ## ---- parse clusters argument into numeric node IDs ----------------------
  if (missing(clusters) || is.null(clusters)) {
    ids <- seq_len(n_nodes)
  } else {
    if (is.list(clusters)) clusters <- unlist(clusters, use.names = FALSE)
    if (is.character(clusters)) {
      ids <- as.integer(sub("^c_", "", clusters))
    } else {
      ids <- as.integer(clusters)
    }
    ids <- ids[is.finite(ids) & ids > 0L & ids <= n_nodes]
    ids <- sort(unique(ids))
  }

  if (length(ids) < 2L) {
    stop("Need at least two valid cluster IDs to compute distances.")
  }

  ## ---- membership matrix (plots × selected clusters) ----------------------
  # membership is Plot.cluster > 0 (works for binary and rel_cover)
  G <- PC[, ids, drop = FALSE]
  G <- Matrix::Matrix(G > 0, sparse = TRUE)
  G <- G * 1

  labels <- paste0("c_", ids)
  colnames(G) <- labels

  ## ---- compute phi similarity matrix via crossproducts --------------------
  # a_ij = number of plots where both clusters present
  A <- as.matrix(Matrix::crossprod(G))   # clusters × clusters
  P <- as.numeric(Matrix::colSums(G))    # membership counts per cluster
  N <- n_plots

  Pi <- matrix(P, nrow = length(P), ncol = length(P))
  Pj <- t(Pi)

  a <- A
  b <- Pi - a
  c <- Pj - a
  d <- N - a - b - c

  den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
  phi <- (a * d - b * c) / den
  phi[!is.finite(phi) | den <= 0] <- 0
  diag(phi) <- 1

  ## ---- distance: always 1 - phi ------------------------------------------
  Dmat <- 1 - phi
  diag(Dmat) <- 0
  dimnames(Dmat) <- list(labels, labels)

  stats::as.dist(Dmat)
}
