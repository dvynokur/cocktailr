#' @keywords internal
.top_clusters_at_phi <- function(Cluster.merged, Cluster.height, phi) {
  CM <- as.matrix(Cluster.merged); H <- as.numeric(Cluster.height)
  idx <- which(H >= phi); if (!length(idx)) return(integer(0))
  children <- unique(as.integer(CM[idx, ][CM[idx, ] > 0]))
  sort(setdiff(idx, intersect(idx, children)))
}

#' Species -> group codes for one cut
#' @param x cocktail object
#' @param phi numeric cut (e.g., 0.3)
#' @param not_assigned label for species not captured
#' @return tibble with species and group code
#' @export
cut_groups <- function(x, phi, not_assigned = "0") {
  CS <- x$Cluster.species; CM <- x$Cluster.merged; H <- x$Cluster.height
  species <- colnames(CS)
  top_ids <- .top_clusters_at_phi(CM, H, phi)
  codes <- rep(not_assigned, ncol(CS)); names(codes) <- species
  if (length(top_ids)) {
    ord <- order(vapply(top_ids, function(i) min(which(CS[i,] == 1)), 1L))
    top_ids <- top_ids[ord]
    for (k in seq_along(top_ids)) {
      s <- which(CS[top_ids[k], ] == 1)
      codes[s] <- sprintf("phi%.2f_G%03d", phi, k)
    }
  }
  tibble::tibble(species = species, group = unname(codes))
}
