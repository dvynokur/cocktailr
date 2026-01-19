#' Distance between Cocktail clusters based on species fidelity (phi)
#'
#' @description
#' Compute a phi-based distance between a set of Cocktail clusters (internal nodes)
#' using species fidelity profiles taken directly from \code{x$Species.cluster.phi}.
#'
#' For each pair of clusters \eqn{A,B}, similarity is computed twice (\eqn{A\to B}
#' and \eqn{B\to A}) and then averaged. In the directional comparison \eqn{A\to B},
#' only species with \eqn{\phi(s,A) \ge \mathrm{phi\_threshold}} are used; the
#' corresponding \eqn{\phi(s,B)} values are taken as-is (can be negative). The
#' final similarity is:
#' \deqn{
#' \mathrm{sim}(A,B)=\frac{1}{2}\left[\mathrm{sim}(A\to B)+\mathrm{sim}(B\to A)\right],
#' }
#' and the distance is \eqn{d(A,B)=1-\mathrm{sim}(A,B)}.
#'
#' Similarities are computed using one of three methods: \code{"dot"},
#' \code{"cosine"}, or \code{"pearson"}.
#'
#' Optionally, a signed power transform is applied to phi values prior to computing
#' similarities to downweight near-zero values:
#' \deqn{\tilde{\phi}=\mathrm{sign}(\phi)\,|\phi|^{\gamma}}
#' with \eqn{\gamma = 2}.
#'
#' @param x A \code{"cocktail"} object (result of \code{\link{cocktail_cluster}}),
#'   containing at least:
#'   \itemize{
#'     \item \code{Cluster.species} - nodes x species (0/1),
#'     \item \code{Cluster.height} - numeric vector of merge phi for each node, and
#'     \item \code{Species.cluster.phi} - species x nodes phi-matrix.
#'   }
#'   Note: \code{Species.cluster.phi} is only present if
#'   \code{species_cluster_phi = TRUE} was used in \code{\link{cocktail_cluster}}.
#'
#' @param clusters Optional cluster identifiers (nodes) to be compared. Can be a
#'   numeric vector of node indices (e.g. \code{c(12, 27)}) or a character vector
#'   of node labels (e.g. \code{c("c_12", "c_27")}). Each element refers to a
#'   single internal node; no grouping/union is performed.
#'   If \code{NULL} or missing, all internal nodes
#'   \code{1:nrow(x$Cluster.species)} are candidates.
#'
#' @param min_phi Numeric scalar; minimum merge phi (cluster height) required for a
#'   node to be retained. Default \code{0.2}. Nodes whose corresponding
#'   \code{x$Cluster.height} value is \emph{strictly less} than \code{min_phi}
#'   are dropped before computing the distance. If fewer than two nodes remain,
#'   an error is raised.
#'
#' @param method Similarity method used to compare cluster fidelity profiles on
#'   the selected species set:
#'   \itemize{
#'     \item \code{"dot"}: mean elementwise product (mean of \code{A * B}). Since
#'       phi values are in \eqn{[-1, 1]}, similarity is also in \eqn{[-1, 1]}.
#'     \item \code{"cosine"}: cosine similarity
#'       \eqn{\frac{A \cdot B}{\|A\|\|B\|}} (returns 0 if a norm is 0).
#'     \item \code{"pearson"}: Pearson correlation (returns 0 if undefined, if
#'       variance is 0, or if fewer than 2 species are available).
#'   }
#'
#' @param phi_threshold Numeric scalar; lower threshold applied to the \emph{source}
#'   cluster in each directional comparison. For \eqn{A\to B}, only species with
#'   \eqn{\phi(s,A) \ge \mathrm{phi\_threshold}} are used. Default \code{-1}.
#'   Practical settings:
#'   \itemize{
#'     \item \code{-1}: include the full profile (includes \eqn{\phi=-1} as well).
#'     \item \code{0}: use only non-negative (positive or zero) fidelities in the source cluster.
#'     \item \code{0.2}: focus on meaningfully associated species in the source cluster.
#'   }
#'   Clusters with no species meeting this threshold are dropped before distance
#'   computation (with a warning).
#'
#' @param drop_nested_clusters Logical. If \code{TRUE}, clusters whose species
#'   sets are strict subsets of other retained clusters are dropped (based
#'   on \code{x$Cluster.species}). If \code{FALSE} (default), nested clusters
#'   are retained.
#'
#' @param power_transform Logical. If \code{TRUE} (default), apply a signed power
#'   transform with \eqn{\gamma=2} to phi values before computing similarities:
#'   \eqn{\tilde{\phi}=\mathrm{sign}(\phi)|\phi|^{2}}. If \code{FALSE}, use raw phi.
#'
#' @details
#' Let \eqn{\phi(s,k)} be the phi-coefficient between species \eqn{s} and cluster
#' \eqn{k} stored in \code{x$Species.cluster.phi}. For clusters \eqn{A,B}, define
#' the directional species subset \eqn{S_A=\{s:\phi(s,A)\ge \mathrm{phi\_threshold}\}}
#' and compute \eqn{\mathrm{sim}(A\to B)} by the chosen \code{method} comparing
#' \eqn{\phi(S_A,A)} to \eqn{\phi(S_A,B)} (optionally after power transformation).
#'
#' Distances are returned as \eqn{d(A,B)=1-\mathrm{sim}(A,B)}. For \code{"cosine"}
#' and \code{"pearson"}, similarities lie in \eqn{[-1, 1]}, hence distances lie
#' in \eqn{[0, 2]}. For \code{"dot"} as defined here (mean elementwise product),
#' similarity also lies in \eqn{[-1, 1]} and distance in \eqn{[0, 2]}.
#'
#' @return A \code{\link[stats]{dist}} object of pairwise distances between clusters.
#'
#' @seealso \code{\link{cocktail_cluster}}, \code{\link{clusters_at_cut}}
#' @importFrom stats as.dist cor sd
#' @export
cluster_phi_dist <- function(
    x,
    clusters = NULL,
    min_phi = 0.2,
    method  = c("dot", "cosine", "pearson"),
    phi_threshold = -1,
    drop_nested_clusters = FALSE,
    power_transform = TRUE
) {
  phi_power <- 2

  ## ---- basic checks -------------------------------------------------------
  if (!is.list(x) || !"Cluster.species" %in% names(x)) {
    stop("`x` must be a Cocktail object with a `Cluster.species` component.")
  }

  if (!"Species.cluster.phi" %in% names(x) || is.null(x$Species.cluster.phi)) {
    stop(
      "`x$Species.cluster.phi` is not available.\n",
      "Recompute the Cocktail clustering with `species_cluster_phi = TRUE`, e.g.:\n",
      "  x <- cocktail_cluster(vegmatrix, species_cluster_phi = TRUE, ...)\n",
      "and then call `cluster_phi_dist(x, ...)`."
    )
  }

  if (!"Cluster.height" %in% names(x) || is.null(x$Cluster.height)) {
    stop("`x$Cluster.height` is not available; cannot apply `min_phi` filtering.")
  }

  if (!is.numeric(min_phi) || length(min_phi) != 1L || is.na(min_phi)) {
    stop("`min_phi` must be a single numeric value.")
  }
  if (!is.numeric(phi_threshold) || length(phi_threshold) != 1L || is.na(phi_threshold)) {
    stop("`phi_threshold` must be a single numeric value.")
  }
  if (!is.logical(drop_nested_clusters) || length(drop_nested_clusters) != 1L || is.na(drop_nested_clusters)) {
    stop("`drop_nested_clusters` must be a single logical value (TRUE/FALSE).")
  }
  if (!is.logical(power_transform) || length(power_transform) != 1L || is.na(power_transform)) {
    stop("`power_transform` must be a single logical value (TRUE/FALSE).")
  }

  method <- match.arg(method)

  CS  <- x$Cluster.species       # nodes x species (0/1)
  Phi <- x$Species.cluster.phi   # species x nodes
  H   <- x$Cluster.height        # length n_nodes

  if (!is.matrix(CS))  stop("`x$Cluster.species` must be a matrix.")
  if (!is.matrix(Phi)) stop("`x$Species.cluster.phi` must be a matrix.")

  n_nodes  <- nrow(CS)
  sp_names <- colnames(CS)
  if (is.null(sp_names)) stop("`x$Cluster.species` must have species column names.")
  if (is.null(rownames(Phi))) stop("`x$Species.cluster.phi` must have species row names.")
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
  Phi <- Phi[sp_names, , drop = FALSE]  # reorder rows to match CS cols
  storage.mode(Phi) <- "double"

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

  if (!length(ids)) {
    stop("No valid cluster IDs found in `clusters` (after filtering to 1..", n_nodes, ").")
  }

  ## ---- filter by min_phi (Cluster.height threshold) -----------------------
  keep_phi <- H[ids] >= min_phi
  if (!any(keep_phi)) {
    stop("No clusters with `Cluster.height >= min_phi` (", min_phi, ") among the requested nodes.")
  }
  if (!all(keep_phi)) {
    dropped <- paste0("c_", ids[!keep_phi])
    warning("Dropping clusters below `min_phi` (", min_phi, "): ", paste(dropped, collapse = ", "))
    ids <- ids[keep_phi]
  }

  ## ---- optionally drop nested clusters (subset species sets) --------------
  if (isTRUE(drop_nested_clusters)) {
    sp_idx_list <- lapply(ids, function(k) which(CS[k, ] > 0))
    sizes <- vapply(sp_idx_list, length, integer(1))
    ord <- order(sizes, decreasing = TRUE)

    kept_pos <- integer(0)
    for (ii in ord) {
      s_i <- sp_idx_list[[ii]]
      if (!length(kept_pos)) {
        kept_pos <- c(kept_pos, ii)
        next
      }
      is_subset <- FALSE
      for (jj in kept_pos) {
        s_j <- sp_idx_list[[jj]]
        if (length(s_j) < length(s_i)) next
        if (length(intersect(s_i, s_j)) == length(s_i)) {
          is_subset <- TRUE
          break
        }
      }
      if (!is_subset) kept_pos <- c(kept_pos, ii)
    }

    kept_pos <- sort(kept_pos)
    ids_kept <- ids[kept_pos]

    if (length(ids_kept) < length(ids)) {
      warning(
        "Dropping nested clusters (subset topological species sets): ",
        paste(setdiff(paste0("c_", ids), paste0("c_", ids_kept)), collapse = ", ")
      )
    }
    ids <- ids_kept
  }

  if (length(ids) < 2L) {
    stop("Fewer than two clusters remain after filtering; cannot build a distance matrix.")
  }

  ## ---- map node IDs to columns of Phi -------------------------------------
  node_colnames <- colnames(Phi)
  if (is.null(node_colnames)) {
    stop("`Species.cluster.phi` must have column names (node IDs, e.g. 'c_1').")
  }

  cols_needed <- paste0("c_", ids)
  col_ids     <- match(cols_needed, node_colnames)
  if (anyNA(col_ids)) {
    missing_cols <- cols_needed[is.na(col_ids)]
    stop("`Species.cluster.phi` is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  ## ---- extract full profiles (cluster x species) ---------------------------
  labels <- paste0("c_", ids)
  P <- t(Phi[, col_ids, drop = FALSE])   # rows = clusters, cols = species
  dimnames(P) <- list(labels, sp_names)

  ## ---- drop clusters with no species meeting threshold --------------------
  eligible <- rowSums(P >= phi_threshold, na.rm = TRUE) > 0L
  if (!any(eligible)) {
    stop("No clusters have any species with phi >= ", phi_threshold,
         " (after `min_phi`/`clusters`/nested filtering).")
  }
  if (!all(eligible)) {
    warning(
      "Dropping clusters with no species having phi >= ", phi_threshold, ": ",
      paste(rownames(P)[!eligible], collapse = ", ")
    )
    P <- P[eligible, , drop = FALSE]
  }

  if (nrow(P) < 2L) {
    stop("Fewer than two clusters remain after `phi_threshold` filtering; cannot build distances.")
  }

  ## ---- similarity helpers -------------------------------------------------
  .transform <- function(v) {
    if (!power_transform) return(v)
    sign(v) * (abs(v) ^ phi_power)
  }

  sim_method <- function(a_raw, b_raw) {
    ok <- is.finite(a_raw) & is.finite(b_raw)
    if (!any(ok)) return(0)

    a <- .transform(a_raw[ok])
    b <- .transform(b_raw[ok])

    if (method == "dot") {
      return(mean(a * b))
    }

    if (method == "cosine") {
      na <- sqrt(sum(a * a))
      nb <- sqrt(sum(b * b))
      if (!is.finite(na) || !is.finite(nb) || na <= 0 || nb <= 0) return(0)
      return(sum(a * b) / (na * nb))
    }

    # pearson
    if (length(a) < 2L) return(0)
    if (stats::sd(a) <= 0 || stats::sd(b) <= 0) return(0)
    cc <- suppressWarnings(stats::cor(a, b, method = "pearson"))
    if (!is.finite(cc) || is.na(cc)) return(0)
    cc
  }

  sim_dir <- function(a_raw, b_raw) {
    # A -> B: keep species where a_raw >= phi_threshold (inclusive)
    idx <- is.finite(a_raw) & is.finite(b_raw) & (a_raw >= phi_threshold)
    if (!any(idx)) return(0)
    sim_method(a_raw[idx], b_raw[idx])
  }

  sim_sym <- function(a_raw, b_raw) {
    0.5 * (sim_dir(a_raw, b_raw) + sim_dir(b_raw, a_raw))
  }

  ## ---- compute distance matrix -------------------------------------------
  n_cl <- nrow(P)
  Dmat <- matrix(0, n_cl, n_cl, dimnames = list(rownames(P), rownames(P)))

  for (i in 1L:(n_cl - 1L)) {
    ai <- P[i, ]
    for (j in (i + 1L):n_cl) {
      bj <- P[j, ]
      s  <- sim_sym(ai, bj)
      Dmat[i, j] <- Dmat[j, i] <- (1 - s)
    }
  }

  stats::as.dist(Dmat)
}
