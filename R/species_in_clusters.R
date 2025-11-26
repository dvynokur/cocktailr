#' Species in Cocktail clusters (topology and optional φ)
#'
#' @description
#' Return diagnostic species sets for selected Cocktail nodes or node unions.
#' Works with a Cocktail object produced by \code{cocktail_cluster()} and uses:
#' \itemize{
#'   \item \code{x$Cluster.species} — topological species per node (binary or
#'         positive weights; any value > 0 is treated as membership),
#'   \item optionally \code{x$Species.cluster.phi} — species–cluster \eqn{\phi}
#'         (if \code{species_cluster_phi = TRUE}).
#' }
#'
#' The \code{labels} argument can be either:
#' \itemize{
#'   \item a vector of node IDs / labels, where each element is treated as a
#'         separate node; or
#'   \item a \strong{list}, where each element is a vector of node IDs / labels,
#'         and each list element is treated as a \strong{union group} of those
#'         nodes (species are the union across the nodes in that group).
#' }
#'
#' Behaviour:
#'
#' \describe{
#'
#' \item{\code{species_cluster_phi = FALSE}}{
#'   For each group, return the union of its \strong{topological species}
#'   (from \code{x$Cluster.species}). Output: named list of character vectors
#'   of species names.
#' }
#'
#' \item{\code{species_cluster_phi = TRUE}}{
#'   Requires \code{x$Species.cluster.phi}. For each group (node or union),
#'   define group-level species–group \eqn{\phi(s,g)} as the \strong{maximum}
#'   across the nodes in that group:
#'   \deqn{\phi(s,g) = \max_{k \in \text{group}} \phi(s,k)}
#'   Negative values are set to 0.
#'
#'   Filtering:
#'   \itemize{
#'     \item If \code{min_phi} is \code{NULL}, the function returns all
#'           \strong{topological species} of the group (from
#'           \code{Cluster.species}), with their group-level \eqn{\phi(s,g)}.
#'     \item If \code{min_phi} is provided (numeric), it returns \emph{all}
#'           species whose group-level \eqn{\phi(s,g) \ge \code{min_phi}},
#'           regardless of whether they are topological or not.
#'   }
#'
#'   Species are sorted by \eqn{\phi} decreasing if \code{sort_desc = TRUE}.
#'   If \code{top_k} is non-NULL, only the top \code{k} species (after sorting
#'   and filtering) are kept.
#'
#'   Output: named list of data frames, one per group, with columns
#'   \code{species} and \code{phi}.
#' }
#'
#' }
#'
#' If \code{species_cluster_phi = TRUE} but \code{x$Species.cluster.phi} is
#' missing, the function issues a warning, suggests recomputing
#' \code{cocktail_cluster(..., species_cluster_phi = TRUE)}, and falls back to
#' returning only topological species (as if \code{species_cluster_phi = FALSE}).
#'
#' @param x A Cocktail object (a list containing at least
#'   \code{Cluster.species}; usually from \code{cocktail_cluster()}).
#' @param labels Selection of clusters. Can be:
#'   \itemize{
#'     \item Character labels like \code{"c_2127"} or integer node IDs
#'           (e.g. \code{2127}), where each element defines one group; or
#'     \item A \strong{list} of such vectors, e.g.
#'           \code{list(10, c(11,13), c(12,15), 14)}, where each element
#'           defines a union group of nodes.
#'   }
#'   Must be supplied.
#' @param species_cluster_phi Logical; if \code{FALSE} (default), return only
#'   topological species per group. If \code{TRUE}, use
#'   \code{x$Species.cluster.phi} to attach \eqn{\phi} values as described.
#' @param min_phi Numeric or \code{NULL}. Used only when
#'   \code{species_cluster_phi = TRUE}:
#'   \itemize{
#'     \item \code{NULL} (default): list only topological species, with their
#'           group-level \eqn{\phi};
#'     \item numeric: list all species with \eqn{\phi(s,g) \ge \code{min_phi}}.
#'   }
#' @param top_k Optional integer ≥ 1. When \code{species_cluster_phi = TRUE},
#'   keep only the top \code{k} species per group after sorting; default
#'   \code{NULL} keeps all.
#' @param sort_desc Logical; when \code{species_cluster_phi = TRUE}, sort
#'   species by \eqn{\phi} decreasing if \code{TRUE} (default).
#'
#' @return
#' \itemize{
#'   \item If \code{species_cluster_phi = FALSE}: a named list of character
#'         vectors (species names per group).
#'   \item If \code{species_cluster_phi = TRUE}: a named list of data frames
#'         (one per group) with columns \code{species} and \code{phi}.
#' }
#'
#' @seealso \code{\link{cocktail_cluster}}, \code{\link{assign_releves}}
#' @export
species_in_clusters <- function(
    x,
    labels,
    species_cluster_phi = FALSE,
    min_phi             = NULL,
    top_k               = NULL,
    sort_desc           = TRUE
) {
  ## ---- helpers ------------------------------------------------------------
  .is_cocktail <- function(obj) {
    is.list(obj) && "Cluster.species" %in% names(obj)
  }

  # Normalise a vector of labels (character or numeric) to integer node IDs
  .norm_label_vector_to_ids <- function(lab) {
    if (is.null(lab)) return(integer(0))
    if (is.character(lab)) {
      out <- sub("^c_", "", lab)
      as.integer(out)
    } else if (is.numeric(lab) || is.integer(lab)) {
      as.integer(lab)
    } else {
      stop("Labels must be character like 'c_123' or integer node IDs.")
    }
  }

  # Normalise labels into a list of node-ID vectors and group names
  .normalize_labels_to_groups <- function(labels) {
    if (is.null(labels)) {
      stop("Please supply `labels` (vector or list) for a Cocktail object.")
    }

    if (is.list(labels)) {
      grp_ids <- lapply(labels, .norm_label_vector_to_ids)
      grp_names <- vapply(
        grp_ids,
        function(ids) paste0("c_", paste(ids, collapse = "_")),
        character(1)
      )
      list(ids = grp_ids, names = grp_names)
    } else {
      ids <- .norm_label_vector_to_ids(labels)
      grp_ids   <- lapply(as.list(ids), identity)
      grp_names <- paste0("c_", ids)
      list(ids = grp_ids, names = grp_names)
    }
  }

  ## ---- checks / setup -----------------------------------------------------
  if (!.is_cocktail(x)) {
    stop("`x` must be a Cocktail object with a `Cluster.species` component.")
  }

  CS <- x$Cluster.species
  sp_names <- colnames(CS)
  if (is.null(sp_names)) {
    stop("`Cluster.species` must have column (species) names.")
  }

  # make sure we can safely compare numerically
  CS_num <- as.matrix(CS)
  storage.mode(CS_num) <- "double"

  n_nodes <- nrow(CS_num)

  labs <- .normalize_labels_to_groups(labels)
  grp_ids   <- labs$ids
  grp_names <- labs$names

  # validate node IDs and drop invalid ones per group
  valid_groups <- logical(length(grp_ids))
  for (i in seq_along(grp_ids)) {
    ids <- grp_ids[[i]]
    bad <- is.na(ids) | ids < 1L | ids > n_nodes
    if (any(bad)) {
      warning(
        "Dropping invalid node IDs in group ", grp_names[i], ": ",
        paste(ids[bad], collapse = ", ")
      )
      ids <- ids[!bad]
      grp_ids[[i]] <- ids
    }
    valid_groups[i] <- length(ids) > 0L
  }

  if (!any(valid_groups)) {
    return(setNames(list(), character(0)))
  }

  grp_ids   <- grp_ids[valid_groups]
  grp_names <- grp_names[valid_groups]

  # If φ requested but matrix missing -> warn and fall back to topology only
  if (isTRUE(species_cluster_phi) &&
      (!("Species.cluster.phi" %in% names(x)) || is.null(x$Species.cluster.phi))) {
    warning(
      "species_cluster_phi = TRUE requested, but x$Species.cluster.phi is not present.\n",
      "Recompute cocktail_cluster(..., species_cluster_phi = TRUE) to enable φ output.\n",
      "Returning only topological species (species_cluster_phi = FALSE)."
    )
    species_cluster_phi <- FALSE
  }

  ## ---- Case 1: topology only ---------------------------------------------
  if (!isTRUE(species_cluster_phi)) {
    out <- vector("list", length(grp_ids))
    for (i in seq_along(grp_ids)) {
      ids <- grp_ids[[i]]
      # treat any value > 0 as membership (works for 0/1 and relative covers)
      sp_logical <- colSums(CS_num[ids, , drop = FALSE] > 0) > 0L
      out[[i]] <- sp_names[sp_logical]
    }
    names(out) <- grp_names
    return(out)
  }

  ## ---- Case 2: φ-based output --------------------------------------------
  Phi <- x$Species.cluster.phi
  if (is.null(Phi) || !is.matrix(Phi)) {
    stop("Internal error: species_cluster_phi = TRUE but `x$Species.cluster.phi` is missing or not a matrix.")
  }
  if (is.null(rownames(Phi))) {
    stop("`x$Species.cluster.phi` must have species row names.")
  }

  # Align species: restrict both CS and Phi to common species
  common <- intersect(sp_names, rownames(Phi))
  if (!length(common)) {
    stop("No overlapping species between `Cluster.species` and `Species.cluster.phi`.")
  }
  CS_num <- CS_num[, common, drop = FALSE]
  Phi    <- Phi[common, , drop = FALSE]
  sp_names <- common

  # normalise min_phi
  if (!is.null(min_phi)) {
    if (!is.numeric(min_phi) || !is.finite(min_phi)) {
      min_phi <- NULL
    } else {
      min_phi <- max(0, min(1, as.numeric(min_phi)))
    }
  }
  do_top_k <- !is.null(top_k) && is.finite(top_k) && top_k > 0

  node_colnames <- colnames(Phi)
  if (is.null(node_colnames)) {
    stop("`Species.cluster.phi` must have column names (node IDs, e.g. 'c_1').")
  }

  .node_to_col_idx <- function(k) {
    cand1 <- paste0("c_", k)
    if (cand1 %in% node_colnames) return(match(cand1, node_colnames))
    cand2 <- as.character(k)
    if (cand2 %in% node_colnames) return(match(cand2, node_colnames))
    NA_integer_
  }

  out <- vector("list", length(grp_ids))

  for (i in seq_along(grp_ids)) {
    ids <- grp_ids[[i]]

    # topological species union for this group (> 0 membership)
    sp_logical <- colSums(CS_num[ids, , drop = FALSE] > 0) > 0L
    sp_topo <- sp_names[sp_logical]

    if (!length(ids)) {
      out[[i]] <- data.frame(
        species = character(0),
        phi     = numeric(0),
        stringsAsFactors = FALSE
      )
      next
    }

    col_idx <- vapply(ids, .node_to_col_idx, integer(1L))
    if (anyNA(col_idx)) {
      missing_ids <- ids[is.na(col_idx)]
      stop("`Species.cluster.phi` is missing columns for nodes: ",
           paste(missing_ids, collapse = ", "))
    }

    Phi_g <- Phi[, col_idx, drop = FALSE]   # species × nodes_in_group
    phi_s <- apply(Phi_g, 1L, max)         # group-level phi per species
    phi_s[phi_s < 0] <- 0

    # Determine which species to keep
    if (is.null(min_phi)) {
      # only topological species, regardless of phi magnitude
      keep_species <- sp_topo
    } else {
      keep_species <- names(phi_s)[phi_s >= min_phi]
    }

    if (!length(keep_species)) {
      out[[i]] <- data.frame(
        species = character(0),
        phi     = numeric(0),
        stringsAsFactors = FALSE
      )
      next
    }

    vals <- phi_s[keep_species]
    df <- data.frame(
      species = keep_species,
      phi     = as.numeric(vals),
      stringsAsFactors = FALSE
    )

    if (isTRUE(sort_desc)) {
      df <- df[order(df$phi, decreasing = TRUE), , drop = FALSE]
    }
    if (do_top_k && nrow(df) > top_k) {
      df <- df[seq_len(top_k), , drop = FALSE]
    }
    rownames(df) <- NULL
    out[[i]] <- df
  }

  names(out) <- grp_names
  out
}
