#' Distance between Cocktail clusters based on species fidelity (φ)
#'
#' @description
#' Compute a φ-based distance between a set of Cocktail clusters (internal nodes).
#' For each cluster, a species fidelity profile is taken directly from
#' \code{x$Species.cluster.phi}, restricted to the cluster's topological species
#' set. Pairwise distances between clusters are then derived from these profiles
#' via a symmetric similarity measure.
#'
#' @param x A \code{"cocktail"} object (result of \code{\link{cocktail_cluster}}),
#'   containing at least:
#'   \itemize{
#'     \item \code{Cluster.species} — nodes × species (0/1),
#'     \item \code{Cluster.height} — numeric vector of merge φ for each node, and
#'     \item \code{Species.cluster.phi} — species × nodes Option A φ-matrix.
#'   }
#'   Note: \code{Species.cluster.phi} is only present if
#'   \code{species_cluster_phi = TRUE} was used in \code{\link{cocktail_cluster}}.
#'
#' @param clusters Optional cluster identifiers (nodes) to be compared. Can be a
#'   numeric vector of node indices (e.g. \code{c(12, 27)}) or a character vector
#'   of node labels (e.g. \code{c("c_12", "c_27")}). Each element refers to a
#'   single internal node; no grouping/union is performed.
#'   If \code{NULL} or missing, **all** internal nodes
#'   \code{1:nrow(x$Cluster.species)} are candidates.
#'
#' @param min_phi Numeric scalar; minimum merge φ (cluster height) required for a
#'   node to be retained. Default \code{0.2}. Nodes whose corresponding
#'   \code{x$Cluster.height} value is \emph{strictly less} than \code{min_phi}
#'   are dropped before computing the φ-based distance. If, after applying
#'   \code{clusters} and \code{min_phi}, fewer than two nodes remain, an error
#'   is raised.
#'
#' @details
#' Let \eqn{\phi(s, k)} be the Option A φ-coefficient between species \eqn{s}
#' and cluster \eqn{k}, stored in \code{x$Species.cluster.phi}. For each
#' cluster \eqn{k}:
#'
#' \enumerate{
#'   \item Its **topological species set** \eqn{S_k} is defined as the set of
#'         species with membership 1 in \code{x$Cluster.species[k,]}.
#'   \item A species fidelity profile \eqn{\phi_k(s)} is taken as
#'         \eqn{\phi(s,k)} for \eqn{s \in S_k}, and 0 for species outside
#'         \eqn{S_k}. All negative φ values are set to 0 beforehand, so only
#'         positive fidelity contributes.
#' }
#'
#' Only nodes with \code{x$Cluster.height >= min_phi} are retained for distance
#' computation (after any user-specified filtering via \code{clusters}).
#'
#' For any pair of clusters \eqn{A,B}, define:
#'
#' \deqn{
#'   \mathrm{sim}(A \to B) =
#'     \frac{\sum_{s \in S_A} \phi_B(s)}
#'          {\sum_{s \in S_A} \phi_A(s)}
#' }
#'
#' with the convention that if the denominator is 0 (no positive φ for
#' \eqn{A}), then \eqn{\mathrm{sim}(A\to B) = 0}. The symmetric similarity is
#'
#' \deqn{
#'   \mathrm{sim}_{\mathrm{sym}}(A,B) =
#'     \frac{1}{2}\left[ \mathrm{sim}(A\to B) + \mathrm{sim}(B\to A) \right],
#' }
#'
#' and the distance is \eqn{d(A,B) = 1 - \mathrm{sim}_{\mathrm{sym}}(A,B)}.
#'
#' @return
#' A \code{\link[stats]{dist}} object of pairwise distances between the clusters.
#' No additional attributes are attached; you can pass the result directly to
#' \code{\link[stats]{hclust}} or other clustering/ordination functions.
#'
#' @seealso \code{\link{cocktail_cluster}}, \code{\link{clusters_at_cut}}
#' @importFrom stats as.dist
#' @export
cluster_phi_dist <- function(
    x,
    clusters = NULL,
    min_phi = 0.2
) {
  ## ---- basic checks -------------------------------------------------------
  if (!is.list(x) || !"Cluster.species" %in% names(x)) {
    stop("`x` must be a Cocktail object with a `Cluster.species` component.")
  }

  # Explicit, user-friendly check for Species.cluster.phi
  if (!"Species.cluster.phi" %in% names(x) || is.null(x$Species.cluster.phi)) {
    stop(
      "`x$Species.cluster.phi` is not available.\n",
      "Recompute the Cocktail clustering with `species_cluster_phi = TRUE`, e.g.:\n",
      "  x <- cocktail_cluster(vegmatrix, species_cluster_phi = TRUE, ...)\n",
      "and then call `cluster_phi_dist(x, ...)`."
    )
  }

  # Need Cluster.height to apply min_phi
  if (!"Cluster.height" %in% names(x) || is.null(x$Cluster.height)) {
    stop("`x$Cluster.height` is not available; cannot apply `min_phi` filtering.")
  }

  if (!is.numeric(min_phi) || length(min_phi) != 1L || is.na(min_phi)) {
    stop("`min_phi` must be a single numeric value.")
  }

  CS   <- x$Cluster.species       # nodes × species (0/1)
  Phi  <- x$Species.cluster.phi   # species × nodes
  H    <- x$Cluster.height        # length n_nodes

  if (!is.matrix(CS))
    stop("`x$Cluster.species` must be a matrix.")
  if (!is.matrix(Phi))
    stop("`x$Species.cluster.phi` must be a matrix.")

  n_nodes  <- nrow(CS)
  sp_names <- colnames(CS)
  if (is.null(sp_names))
    stop("`x$Cluster.species` must have species column names.")

  if (is.null(rownames(Phi)))
    stop("`x$Species.cluster.phi` must have species row names.")

  if (length(H) < n_nodes) {
    stop("`x$Cluster.height` must be a numeric vector of length nrow(x$Cluster.species).")
  }

  ## ---- align species between CS and Phi -----------------------------------
  if (!all(sp_names %in% rownames(Phi))) {
    missing_sp <- setdiff(sp_names, rownames(Phi))
    stop("`Species.cluster.phi` is missing species: ",
         paste(head(missing_sp, 10), collapse = ", "),
         if (length(missing_sp) > 10) " ..." else "")
  }
  Phi <- Phi[sp_names, , drop = FALSE]  # reorder rows

  # keep only non-negative fidelity
  Phi_pos <- pmax(Phi, 0)
  storage.mode(Phi_pos) <- "double"

  ## ---- parse clusters argument into numeric node IDs ----------------------
  if (missing(clusters) || is.null(clusters)) {
    # use all internal nodes
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

  if (!length(ids)) {
    stop("No valid cluster IDs found in `clusters` (after filtering to 1..",
         n_nodes, ").")
  }

  ## ---- filter by min_phi (Cluster.height threshold) -----------------------
  keep_phi <- H[ids] >= min_phi
  if (!any(keep_phi)) {
    stop("No clusters with `Cluster.height >= min_phi` (", min_phi,
         ") among the requested nodes.")
  }
  if (!all(keep_phi)) {
    dropped <- paste0("c_", ids[!keep_phi])
    warning("Dropping clusters below `min_phi` (", min_phi, "): ",
            paste(dropped, collapse = ", "))
    ids <- ids[keep_phi]
  }

  ## ---- map node IDs to columns of Phi_pos ---------------------------------
  node_colnames <- colnames(Phi_pos)
  if (is.null(node_colnames)) {
    stop("`Species.cluster.phi` must have column names (node IDs, e.g. 'c_1').")
  }

  cols_needed <- paste0("c_", ids)
  col_ids     <- match(cols_needed, node_colnames)
  if (anyNA(col_ids)) {
    missing_cols <- cols_needed[is.na(col_ids)]
    stop("`Species.cluster.phi` is missing columns: ",
         paste(missing_cols, collapse = ", "))
  }

  ## ---- build cluster × species φ-matrix (species fidelity profiles) -------
  n_cl   <- length(ids)
  n_sp   <- length(sp_names)
  labels <- paste0("c_", ids)

  M <- matrix(0, nrow = n_cl, ncol = n_sp,
              dimnames = list(labels, sp_names))

  for (i in seq_len(n_cl)) {
    k <- ids[i]

    # topological species set: species with membership 1 in Cluster.species[k,]
    sp_idx <- which(CS[k, ] > 0)

    if (!length(sp_idx)) {
      next  # this cluster will become all zeros; may be dropped later
    }

    col_k <- col_ids[i]
    phi_k <- Phi_pos[sp_idx, col_k]

    M[i, sp_idx] <- phi_k
  }

  ## ---- drop clusters with completely zero profiles (no positive φ) -------
  zero_rows <- which(rowSums(M) <= 0)
  if (length(zero_rows)) {
    if (length(zero_rows) == nrow(M)) {
      stop("All requested clusters have zero positive φ in their ",
           "topological species sets; cannot compute distances.")
    }
    warning("Some clusters have zero positive φ in their topological species ",
            "sets and will be dropped: ",
            paste(rownames(M)[zero_rows], collapse = ", "))
    keep <- setdiff(seq_len(nrow(M)), zero_rows)
    M    <- M[keep, , drop = FALSE]
    labels <- labels[keep]
    ids    <- ids[keep]
  }

  n_cl <- nrow(M)
  if (n_cl < 2L) {
    stop("Fewer than two clusters remain after filtering; cannot build a ",
         "distance matrix.")
  }

  ## ---- compute symmetric φ-based distance matrix --------------------------
  Dmat <- matrix(0, n_cl, n_cl,
                 dimnames = list(labels, labels))

  for (i in 1L:(n_cl - 1L)) {
    Si      <- which(M[i, ] > 0)
    denom_i <- sum(M[i, Si])

    for (j in (i + 1L):n_cl) {
      Sj      <- which(M[j, ] > 0)
      denom_j <- sum(M[j, Sj])

      sim_i_to_j <- if (length(Si) && denom_i > 0) sum(M[j, Si]) / denom_i else 0
      sim_j_to_i <- if (length(Sj) && denom_j > 0) sum(M[i, Sj]) / denom_j else 0

      sim_sym <- 0.5 * (sim_i_to_j + sim_j_to_i)
      d_ij    <- 1 - sim_sym

      Dmat[i, j] <- Dmat[j, i] <- d_ij
    }
  }

  stats::as.dist(Dmat)
}
