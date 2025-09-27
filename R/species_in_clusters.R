#' Species in Cocktail clusters or fuzzy groups
#'
#' @description
#' Return the diagnostic species set for selected Cocktail nodes.
#' Works with either a Cocktail object (membership is binary) or a
#' species-by-nodes \eqn{\phi} matrix produced by \code{cocktail_fuzzy()}.
#'
#' @param x Either:
#'   \itemize{
#'     \item a Cocktail object (a list containing \code{Cluster.species}); or
#'     \item a numeric matrix returned by \code{cocktail_fuzzy()} with species in
#'           rownames and node labels like \code{"c_2127"} in colnames.
#'   }
#' @param labels Character cluster labels such as \code{"c_2127"} or integer node IDs.
#'   If \code{x} is a \eqn{\phi} matrix and \code{labels} is omitted, all columns are used.
#'   If \code{x} is a Cocktail object, \code{labels} must be supplied.
#' @param min_phi Numeric threshold used only when \code{x} is a \eqn{\phi} matrix:
#'   keep species with \eqn{\phi \ge} \code{min_phi}. Default is \code{0.20}.
#' @param top_k Optional integer used only when \code{x} is a \eqn{\phi} matrix:
#'   after filtering by \code{min_phi}, keep the top \code{k} species by \eqn{\phi}
#'   for each requested label. Default \code{NULL} keeps all that pass \code{min_phi}.
#' @param sort_desc Logical; when \code{x} is a \eqn{\phi} matrix, sort species by
#'   \eqn{\phi} decreasing (default \code{TRUE}).
#'
#' @return
#'   \itemize{
#'     \item If \code{x} is a Cocktail object: a named list of character vectors
#'           (species names) for the requested labels.
#'     \item If \code{x} is a \eqn{\phi} matrix: a named list of data frames
#'           (one per label) with columns \code{species} and \code{phi}.
#'   }
#'
#' @seealso \code{\link{cocktail_fuzzy}}
#' @export

species_in_clusters <- function(x, labels = NULL, min_phi = 0.20, top_k = NULL, sort_desc = TRUE) {
  # Helper: detect Cocktail object
  .is_cocktail <- function(obj) {
    is.list(obj) && "Cluster.species" %in% names(obj)
  }

  # Normalize labels
  .norm_labels <- function(lab) {
    if (is.null(lab)) return(NULL)
    if (is.character(lab)) {
      out <- sub("^c_", "", lab)
      as.integer(out)
    } else if (is.numeric(lab) || is.integer(lab)) {
      as.integer(lab)
    } else {
      stop("`labels` must be character like 'c_123' or integer node IDs.")
    }
  }

  # --- Case A: Cocktail object -> return membership species ----------------
  if (.is_cocktail(x)) {
    CS <- x$Cluster.species
    sp_names <- colnames(CS)
    if (is.null(sp_names)) stop("Cluster.species must have column (species) names.")

    if (is.null(labels)) {
      stop("When `x` is a Cocktail object, please supply `labels` (e.g., 'c_2127').")
    }

    node_ids <- .norm_labels(labels)
    n_nodes <- nrow(CS)
    bad <- is.na(node_ids) | node_ids < 1L | node_ids > n_nodes
    if (any(bad)) {
      warning("Dropping invalid labels: ", paste(labels[bad], collapse = ", "))
      labels   <- labels[!bad]
      node_ids <- node_ids[!bad]
    }
    if (!length(node_ids)) return(setNames(list(), character(0)))

    out <- vector("list", length(node_ids))
    for (i in seq_along(node_ids)) {
      nid <- node_ids[i]
      out[[i]] <- sp_names[which(CS[nid, ] == 1L)]
    }
    names(out) <- ifelse(grepl("^c_", labels), labels, paste0("c_", node_ids))
    return(out)
  }

  # --- Case B: φ matrix -> return species + φ data frames ------------------
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be either a Cocktail object or a numeric matrix from cocktail_fuzzy().")
  }
  Phi <- x
  if (is.null(rownames(Phi))) stop("φ matrix must have species in rownames.")
  if (is.null(colnames(Phi))) stop("φ matrix must have node labels in colnames (e.g., 'c_2127').")

  # If labels not supplied, use all columns
  if (is.null(labels)) {
    labels_use <- colnames(Phi)
  } else if (is.character(labels)) {
    lab_char <- ifelse(grepl("^c_\\d+$", labels), labels, paste0("c_", sub("^c_", "", labels)))
    labels_use <- intersect(lab_char, colnames(Phi))
    missing <- setdiff(ifelse(grepl("^c_", labels), labels, paste0("c_", labels)), labels_use)
    if (length(missing)) warning("Dropping unknown labels in φ matrix: ", paste(missing, collapse = ", "))
  } else if (is.numeric(labels) || is.integer(labels)) {
    lab_char <- paste0("c_", as.integer(labels))
    labels_use <- intersect(lab_char, colnames(Phi))
    missing <- setdiff(lab_char, labels_use)
    if (length(missing)) warning("Dropping unknown labels in φ matrix: ", paste(missing, collapse = ", "))
  } else {
    stop("`labels` must be character like 'c_123' or integer node IDs.")
  }

  if (!length(labels_use)) return(setNames(list(), character(0)))

  # φ filters (φ-matrix case only)
  do_min_phi <- is.finite(min_phi)
  do_top_k   <- !is.null(top_k) && is.finite(top_k) && top_k > 0

  out <- vector("list", length(labels_use))
  for (j in seq_along(labels_use)) {
    node <- labels_use[j]
    v <- Phi[, node]
    df <- data.frame(
      species = rownames(Phi),
      phi     = as.numeric(v),
      stringsAsFactors = FALSE
    )
    if (isTRUE(sort_desc)) df <- df[order(df$phi, decreasing = TRUE), , drop = FALSE]
    if (do_min_phi)        df <- df[df$phi >= min_phi, , drop = FALSE]
    if (do_top_k && nrow(df) > top_k) df <- df[seq_len(top_k), , drop = FALSE]
    rownames(df) <- NULL
    out[[j]] <- df
  }
  names(out) <- labels_use
  out
}
