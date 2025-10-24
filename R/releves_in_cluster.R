#' Relevés (plots) passing the m-threshold of Cocktail clusters
#'
#' @description
#' For each requested cluster (node), return the plots that contain at least
#' the node-specific threshold \code{m} species from that cluster’s diagnostic set.
#' The \code{m} threshold is taken directly from the Cocktail object
#' (\code{x$Cluster.info[, 2]}), as determined by the clustering algorithm.
#'
#' @param x A Cocktail object (result of \code{cocktail_cluster()}), containing
#'   at least \code{Cluster.species}, \code{Cluster.info}, and \code{plots}.
#' @param vegmatrix Matrix or data frame with plots in rows and species in columns.
#'   Values \eqn{>} 0 are treated as presences; \code{NA} or \eqn{\le} 0 as absences.
#' @param clusters Character or integer vector giving the cluster IDs to inspect.
#'   Accepts forms like \code{5}, \code{"5"}, or \code{"c_5"}. One or many may be supplied.
#'
#' @return
#' A named list whose names are the **numeric node IDs**.
#' Each element is a character vector of plot names that meet or exceed
#' the stored \code{m} threshold (i.e., contain at least \code{m} member species
#' of that cluster). If no plots qualify, the element is \code{character(0)}.
#'
#' @seealso
#' \code{\link{cocktail_cluster}}, \code{\link{species_in_clusters}}, \code{\link{assign_releves}}
#'
#' @examples
#' vm <- matrix(c(1,0,1,
#'                0,1,0,
#'                1,1,0),
#'              nrow = 3, byrow = TRUE,
#'              dimnames = list(paste0("plot", 1:3),
#'                              c("sp1","sp2","sp3")))
#' res <- cocktail_cluster(vm, progress = FALSE)
#' releves_in_cluster(res, vm, clusters = c("c_1", "c_2"))
#'
#' @export

releves_in_cluster <- function(x, vegmatrix, clusters) {
  ## ---- helpers ----
  .is_cocktail <- function(obj) {
    is.list(obj) && all(c("Cluster.species","Cluster.info") %in% names(obj))
  }
  .norm_clusters <- function(v) {
    if (is.character(v)) as.integer(sub("^c_", "", v))
    else if (is.numeric(v) || is.integer(v)) as.integer(v)
    else stop("`clusters` must be character like 'c_123' or integer node IDs.")
  }

  ## ---- checks ----
  if (!.is_cocktail(x))
    stop("`x` must be a Cocktail object (result of cocktail_cluster()).")

  CS <- x$Cluster.species
  CI <- x$Cluster.info
  if (is.null(CI) || ncol(CI) < 2)
    stop("`x$Cluster.info` must exist with at least two columns (k and m).")

  node_ids <- .norm_clusters(clusters)
  n_nodes <- nrow(CS)
  bad <- is.na(node_ids) | node_ids < 1L | node_ids > n_nodes
  if (any(bad)) {
    warning("Dropping invalid cluster ids: ", paste(clusters[bad], collapse = ", "))
    node_ids <- node_ids[!bad]
  }
  if (!length(node_ids)) return(setNames(list(), character(0)))

  if (!is.matrix(vegmatrix) && !is.data.frame(vegmatrix))
    stop("`vegmatrix` must be a matrix or data.frame (plots × species).")

  vm <- as.matrix(vegmatrix)
  vm[is.na(vm)] <- 0
  vm <- (vm > 0) * 1L
  plots <- rownames(vm)
  if (is.null(plots)) plots <- paste0("plot_", seq_len(nrow(vm)))

  sp_names <- colnames(CS)
  if (is.null(sp_names)) stop("`x$Cluster.species` must have species column names.")
  common <- intersect(colnames(vm), sp_names)
  if (!length(common)) stop("No overlapping species between vegmatrix and Cocktail object.")
  vm <- vm[, common, drop = FALSE]
  CS <- CS[, common, drop = FALSE]

  ## ---- main computation ----
  m_per_node <- as.integer(CI[, 2])
  names(m_per_node) <- as.character(seq_len(nrow(CI)))

  out <- vector("list", length(node_ids))
  for (i in seq_along(node_ids)) {
    nid <- node_ids[i]
    sp_in_cluster <- sp_names[CS[nid, ] == 1L]
    sp_in_cluster <- intersect(sp_in_cluster, colnames(vm))
    if (!length(sp_in_cluster)) {
      out[[i]] <- character(0)
      next
    }

    m_thr <- m_per_node[nid]
    if (!is.finite(m_thr) || m_thr <= 0) m_thr <- 1L

    counts <- rowSums(vm[, sp_in_cluster, drop = FALSE])
    out[[i]] <- plots[counts >= m_thr]
  }

  names(out) <- as.character(node_ids)
  out
}
