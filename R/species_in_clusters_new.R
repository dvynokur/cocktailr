#' Species in Cocktail clusters (with optional unions)
#'
#' @description
#' Return the topological species sets for selected Cocktail nodes.
#' Works only with a Cocktail object and uses its \code{Cluster.species}
#' component (binary membership).
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
#' @param x A Cocktail object (a list containing \code{Cluster.species}).
#' @param labels Selection of clusters. Can be:
#'   \itemize{
#'     \item Character labels like \code{"c_2127"} or integer node IDs
#'           (e.g. \code{2127}), in which case each label defines one group; or
#'     \item A \strong{list} of such vectors, e.g.
#'           \code{list(10, c(11,13), c(12,15), 14)}, where each element
#'           defines a union group of nodes.
#'   }
#'   When \code{x} is a Cocktail object, \code{labels} must be supplied.
#'
#' @return
#' A named list of character vectors. Each element corresponds to one
#' requested node or union group and contains the species names that belong
#' to that group (topological union across the selected nodes).
#'
#' @seealso \code{\link{cocktail_cluster}}
#' @export
species_in_clusters_new <- function(x, labels) {
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
  # For Cocktail: node IDs are integer indices (rows of Cluster.species)
  .normalize_labels_to_groups <- function(labels) {
    if (is.null(labels)) {
      stop("Please supply `labels` (vector or list) for a Cocktail object.")
    }

    if (is.list(labels)) {
      # list of vectors → each element is a group
      grp_ids <- lapply(labels, .norm_label_vector_to_ids)
      grp_names <- vapply(
        grp_ids,
        function(ids) paste0("c_", paste(ids, collapse = "_")),
        character(1)
      )
      list(ids = grp_ids, names = grp_names)
    } else {
      # vector of labels → one group per label
      ids <- .norm_label_vector_to_ids(labels)
      grp_ids <- lapply(as.list(ids), identity)
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

  n_nodes <- nrow(CS)

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

  ## ---- union species per group -------------------------------------------
  out <- vector("list", length(grp_ids))
  for (i in seq_along(grp_ids)) {
    ids <- grp_ids[[i]]
    # union of species across these nodes
    sp_logical <- colSums(CS[ids, , drop = FALSE] == 1L) > 0L
    out[[i]] <- sp_names[sp_logical]
  }
  names(out) <- grp_names
  out
}
