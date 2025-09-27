#' Parent cluster labels at a phi cut
#'
#' @description
#' Return the **topmost** Cocktail nodes whose height is at least the given
#' \eqn{\phi} cut. “Topmost” means nodes at or above the cut that do **not**
#' have an ancestor also at or above the cut.
#'
#' @param x A Cocktail object (list with `Cluster.merged` and `Cluster.height`);
#'   alternatively, a numeric matrix produced by `cocktail_fuzzy()` **that carries**
#'   these two components as attributes (i.e., `attr(x, "Cluster.merged")` and
#'   `attr(x, "Cluster.height")`).
#' @param phi Numeric scalar between 0 and 1 specifying the \eqn{\phi} cut level.
#' @param as_labels Logical; if `TRUE` (default) return character labels like
#'   `"c_2127"`. If `FALSE`, return integer node IDs.
#'
#' @return A vector of cluster identifiers (character when `as_labels = TRUE`,
#'   otherwise integers), sorted increasingly by node ID. May be length zero if
#'   no parent nodes exist at the requested cut.
#'
#' @details
#' This function finds all nodes with `Cluster.height >= phi`, removes any that
#' are descendants of another node that also meets the cut, and returns the
#' remaining parent nodes. When a \eqn{\phi} matrix from `cocktail_fuzzy()` is
#' supplied, the function reads the tree from its attributes.
#'
#' @seealso `assign_releves()`, `species_in_clusters()`, `cocktail_fuzzy()`
#' @export

clusters_at_cut <- function(x, phi, as_labels = TRUE) {
  # Extract tree parts from either Cocktail object or fuzzy matrix with attributes
  if (is.list(x) && all(c("Cluster.merged", "Cluster.height") %in% names(x))) {
    CM <- x$Cluster.merged
    H  <- x$Cluster.height
  } else if (is.matrix(x) &&
             !is.null(attr(x, "Cluster.merged")) &&
             !is.null(attr(x, "Cluster.height"))) {
    CM <- attr(x, "Cluster.merged")
    H  <- attr(x, "Cluster.height")
  } else {
    stop("`x` must be a Cocktail object or a fuzzy matrix with ",
         "attributes 'Cluster.merged' and 'Cluster.height'.")
  }

  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi) || phi < 0 || phi > 1) {
    stop("`phi` must be a single number between 0 and 1.")
  }

  # Nodes at/above the cut
  idx <- which(H >= phi)
  if (!length(idx)) {
    return(if (as_labels) character(0) else integer(0))
  }

  # Children of nodes at/above the cut
  kids_mat <- CM[idx, , drop = FALSE]
  children <- unique(as.integer(kids_mat[kids_mat > 0]))

  # Keep only those nodes in idx that are not a child of another node in idx
  top_idx <- sort(setdiff(idx, intersect(idx, children)))

  if (as_labels) paste0("c_", top_idx) else top_idx
}
