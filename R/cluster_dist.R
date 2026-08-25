#' Distance between Cocktail clusters
#'
#' @description
#' Computes pairwise distances between Cocktail clusters based on binary
#' plot-membership vectors derived from \code{x$Plot.cluster}.
#'
#' Cluster membership is defined as:
#' \code{membership = (Plot.cluster > 0)}. This works for both
#' \code{plot_values = "binary"} and \code{plot_values = "rel_cover"}, because
#' in both cases the stored values are positive for plots satisfying the
#' Cocktail membership rule of a cluster.
#'
#' Two distance methods are available:
#' \itemize{
#'   \item \code{method = "containment"} calculates distance as
#'   \eqn{1 - \max(C_{A|B}, C_{B|A})}, where
#'   \eqn{C_{A|B}} is the proportion of plots satisfying the membership rule
#'   of cluster \eqn{A} that also satisfy the membership rule of cluster
#'   \eqn{B}. This method emphasizes
#'   nested or partly nested plot-membership relationships among clusters.
#'   \item \code{method = "phi"} calculates distance as \eqn{1 - \phi},
#'   where \eqn{\phi} is the phi coefficient between binary plot-membership
#'   vectors. This method emphasizes symmetric co-membership association.
#' }
#'
#' @param x A \code{"cocktail"} object from \code{\link{cocktail_cluster}},
#'   containing at least \code{Plot.cluster} and \code{Cluster.species}.
#'
#' @param clusters Optional cluster identifiers to be compared. Can be numeric
#'   cluster IDs, e.g. \code{c(12, 27)}, or character labels, e.g.
#'   \code{c("c_12", "c_27")}. If \code{NULL} or missing, all clusters
#'   \code{1:nrow(x$Cluster.species)} are used.
#'
#' @param method Character. Distance method to use:
#'   \itemize{
#'     \item \code{"containment"}: distance based on plot-membership
#'     containment, calculated as \eqn{1 - \max(C_{A|B}, C_{B|A})};
#'     \item \code{"phi"}: distance based on the phi coefficient between
#'     binary plot-membership vectors, calculated as \eqn{1 - \phi}.
#'   }
#'
#' @param return Character. Type of output:
#'   \itemize{
#'     \item \code{"dist"}: a \code{\link[stats]{dist}} object;
#'     \item \code{"matrix"}: a symmetric distance matrix;
#'     \item \code{"table"}: a data frame with pairwise cluster relationships
#'     and distances.
#'   }
#'
#' @return Depending on \code{return}, a \code{\link[stats]{dist}} object,
#'   a symmetric matrix, or a data frame with pairwise cluster distances.
#'
#' @details
#' For two clusters \eqn{A} and \eqn{B}, let \eqn{n_A} and \eqn{n_B} be the
#' numbers of plots satisfying their membership rules, and let
#' \eqn{n_{AB}} be the number of plots satisfying both rules.
#'
#' For \code{method = "containment"}, directional containment is calculated as:
#' \deqn{
#' C_{A|B} = \frac{n_{AB}}{n_A}, \qquad
#' C_{B|A} = \frac{n_{AB}}{n_B}.
#' }
#' The similarity is the larger of the two directional containment values:
#' \deqn{
#' \mathrm{sim}(A,B) = \max(C_{A|B}, C_{B|A}),
#' }
#' and the distance is:
#' \deqn{
#' d(A,B) = 1 - \mathrm{sim}(A,B).
#' }
#'
#' For \code{method = "phi"}, the phi coefficient is computed from the 2 x 2
#' table of plot co-membership:
#' \itemize{
#'   \item \eqn{a}: plots where \eqn{A=1} and \eqn{B=1},
#'   \item \eqn{b}: plots where \eqn{A=1} and \eqn{B=0},
#'   \item \eqn{c}: plots where \eqn{A=0} and \eqn{B=1},
#'   \item \eqn{d}: plots where \eqn{A=0} and \eqn{B=0}.
#' }
#' \deqn{
#' \phi = \frac{ad - bc}{\sqrt{(a+c)(b+d)(a+b)(c+d)}}.
#' }
#' Undefined cases with zero denominator are set to \eqn{\phi = 0}.
#'
#' @importFrom stats as.dist
#' @import Matrix
#' @export

cluster_dist <- function(
    x,
    clusters = NULL,
    method = c("containment", "phi"),
    return = c("dist", "matrix", "table")
) {

  ## ---- arguments ----------------------------------------------------------
  method <- match.arg(method)
  return <- match.arg(return)

  ## ---- basic checks -------------------------------------------------------
  if (!is.list(x) || !"Cluster.species" %in% names(x)) {
    stop("`x` must be a Cocktail object with a `Cluster.species` component.")
  }
  if (!"Plot.cluster" %in% names(x) || is.null(x$Plot.cluster)) {
    stop("`x$Plot.cluster` is missing; cannot compute cluster distances.")
  }

  CS <- x$Cluster.species
  PC <- x$Plot.cluster

  if (!is.matrix(CS)) {
    stop("`x$Cluster.species` must be a matrix.")
  }

  # allow base matrices but prefer sparse Matrix objects
  if (!inherits(PC, "Matrix")) {
    PC <- Matrix::Matrix(as.matrix(PC), sparse = TRUE)
  }

  n_nodes <- nrow(CS)
  n_plots <- nrow(PC)

  if (n_plots < 1L) {
    stop("`x$Plot.cluster` has zero plots (rows).")
  }
  if (ncol(PC) < n_nodes) {
    stop("`x$Plot.cluster` must have at least nrow(x$Cluster.species) columns (one per node).")
  }

  ## ---- parse clusters argument into numeric node IDs ----------------------
  if (missing(clusters) || is.null(clusters)) {
    ids <- seq_len(n_nodes)
  } else {
    if (is.list(clusters)) {
      clusters <- unlist(clusters, use.names = FALSE)
    }

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

  ## ---- membership matrix (plots x selected clusters) ----------------------
  # membership is Plot.cluster > 0 (works for binary and rel_cover)
  G <- PC[, ids, drop = FALSE]
  G <- Matrix::Matrix(G > 0, sparse = TRUE)
  G <- G * 1

  labels <- paste0("c_", ids)
  colnames(G) <- labels

  ## ---- shared membership counts ------------------------------------------
  # shared_n_ij = number of plots where both clusters are present
  shared_n <- as.matrix(Matrix::crossprod(G))  # clusters x clusters
  n <- as.numeric(Matrix::colSums(G))          # membership counts per cluster
  N <- n_plots

  ni <- matrix(n, nrow = length(n), ncol = length(n))
  nj <- t(ni)

  ## ---- containment similarity --------------------------------------------
  containment_1_in_2 <- shared_n / ni
  containment_2_in_1 <- shared_n / nj

  containment_1_in_2[!is.finite(containment_1_in_2)] <- 0
  containment_2_in_1[!is.finite(containment_2_in_1)] <- 0

  containment_max <- pmax(containment_1_in_2, containment_2_in_1)
  containment_min <- pmin(containment_1_in_2, containment_2_in_1)

  union_n <- ni + nj - shared_n
  jaccard <- shared_n / union_n
  jaccard[!is.finite(jaccard)] <- 0

  diag(containment_1_in_2) <- 1
  diag(containment_2_in_1) <- 1
  diag(containment_max) <- 1
  diag(containment_min) <- 1
  diag(jaccard) <- 1

  ## ---- phi similarity -----------------------------------------------------
  a <- shared_n
  b <- ni - a
  c <- nj - a
  d <- N - a - b - c

  den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
  phi <- (a * d - b * c) / den
  phi[!is.finite(phi) | den <= 0] <- 0
  diag(phi) <- 1

  ## ---- distance matrix ----------------------------------------------------
  if (method == "containment") {
    Dmat <- 1 - containment_max
  } else {
    Dmat <- 1 - phi
  }

  diag(Dmat) <- 0
  dimnames(Dmat) <- list(labels, labels)

  ## ---- return -------------------------------------------------------------
  if (return == "dist") {
    return(stats::as.dist(Dmat))
  }

  if (return == "matrix") {
    return(Dmat)
  }

  if (return == "table") {
    idx <- which(upper.tri(Dmat), arr.ind = TRUE)

    out <- data.frame(
      cluster_1 = labels[idx[, 1]],
      cluster_2 = labels[idx[, 2]],
      id_1 = ids[idx[, 1]],
      id_2 = ids[idx[, 2]],
      n_1 = n[idx[, 1]],
      n_2 = n[idx[, 2]],
      shared_n = shared_n[idx],
      containment_1_in_2 = containment_1_in_2[idx],
      containment_2_in_1 = containment_2_in_1[idx],
      containment_max = containment_max[idx],
      containment_min = containment_min[idx],
      jaccard = jaccard[idx],
      phi = phi[idx],
      method = method,
      distance = Dmat[idx],
      stringsAsFactors = FALSE
    )

    rownames(out) <- NULL
    return(out)
  }
}
