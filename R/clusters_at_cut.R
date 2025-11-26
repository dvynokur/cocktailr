#' Parent Cocktail clusters at a phi cut
#'
#' @description
#' Return the Cocktail nodes that form the **cluster partition** at a given
#' \eqn{\phi} cut. These are the nodes with height \eqn{\ge \phi} that are
#' **not contained inside any other node** at or above the same cut.
#'
#' In other words, if you slice the Cocktail tree horizontally at height
#' \code{phi}, the function returns exactly the nodes whose branches
#' intersect that cut line.
#'
#' @param x A Cocktail object (list) with components
#'   \code{Cluster.merged} and \code{Cluster.height} as produced by
#'   \code{\link{cocktail_cluster}()}.
#' @param phi Numeric scalar between 0 and 1 specifying the \eqn{\phi} cut level.
#' @param as_labels Logical; if \code{TRUE} (default) return character labels
#'   like \code{"c_12"}. If \code{FALSE}, return integer node IDs.
#'
#' @return
#' A vector of cluster identifiers (character when \code{as_labels = TRUE},
#' otherwise integers), sorted increasingly by node ID. May be length zero if
#' no nodes exist at the requested cut.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Finds all nodes whose height satisfies \code{Cluster.height >= phi}.
#'   \item Looks at their children in \code{Cluster.merged}.
#'   \item Removes any node that appears as a child of another node that also
#'         meets the cut.
#' }
#'
#' The remaining nodes are those that “touch” the cut line without being
#' strictly inside another selected node. They correspond to the natural set
#' of groups at that \eqn{\phi} level and are the same nodes used by functions
#' like \code{\link{assign_releves}()} when only \code{phi_cut} is specified.
#'
#' @seealso \code{\link{assign_releves}}, \code{\link{cocktail_cluster}}
#' @export
clusters_at_cut <- function(x, phi, as_labels = TRUE) {
  if (!is.list(x) || !all(c("Cluster.merged", "Cluster.height") %in% names(x))) {
    stop("`x` must be a Cocktail object with components ",
         "`Cluster.merged` and `Cluster.height` (from cocktail_cluster()).")
  }

  CM <- x$Cluster.merged
  H  <- x$Cluster.height

  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi) || phi < 0 || phi > 1) {
    stop("`phi` must be a single number between 0 and 1.")
  }

  idx <- which(H >= phi)
  if (!length(idx)) {
    return(if (as_labels) character(0) else integer(0))
  }

  kids_mat <- CM[idx, , drop = FALSE]
  children <- unique(as.integer(kids_mat[kids_mat > 0]))

  top_idx <- sort(setdiff(idx, intersect(idx, children)))

  if (as_labels) paste0("c_", top_idx) else top_idx
}
