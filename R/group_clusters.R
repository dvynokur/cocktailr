#' Partition Cocktail nodes (clusters) into k groups
#'
#' @description
#' Partitions **nodes of the Cocktail tree** (i.e., internal clusters) into
#' \code{k_partition} groups using either k-means or hierarchical clustering.
#'
#' Steps:
#' 1) **Select nodes:** keep **all** nodes whose height \eqn{\phi} is
#'    \code{>= phi_cut}. (No “topmost” filtering: ancestors and descendants
#'    that pass the cut are all kept.)
#' 2) **Feature matrix (nodes × plots):** take \code{x$Plot.cluster} columns for
#'    the kept nodes and transpose to get a binary matrix with **rows = nodes**,
#'    **columns = plots**.
#' 3) **Deduplicate identical nodes:** if two nodes have identical plot
#'    memberships, drop duplicates **keeping the rightmost node** (highest column
#'    index in the original \code{Plot.cluster}).
#' 4) **(Optional) feature scaling:** standardize plot-columns when
#'    \code{feature_scale = TRUE} so plots contribute comparably.
#' 5) **Partition nodes:**
#'    - \code{method = "kmeans"}: k-means on the node × plot matrix.
#'    - \code{method = "hclust"}: distances (default Manhattan) + linkage
#'      (default UPGMA/"average") → cut into \code{k_partition} groups.
#'      Group IDs are **relabeled** in dendrogram left-to-right order for readability.
#'
#' @param x Cocktail object with \code{Cluster.height} and \code{Plot.cluster}.
#' @param phi_cut Numeric in (0,1). Keep **all** nodes with height \eqn{\ge} \code{phi_cut}.
#' @param k_partition Positive integer, number of node groups to form.
#' @param method \code{"kmeans"} or \code{"hclust"}.
#' @param feature_scale Logical; if \code{TRUE} (default), standardize plot-columns.
#' @param kmeans_seed Optional integer seed for reproducible k-means.
#' @param dist_method Distance for \code{method="hclust"} (default \code{"manhattan"}).
#' @param clust_method Linkage for \code{method="hclust"} (default \code{"average"}).
#' @param use_parallelDist Logical; if \code{TRUE}, try \code{parallelDist::parDist()};
#'   fallback to \code{stats::dist()}.
#' @param use_fastcluster Logical; if \code{TRUE}, try \code{fastcluster::hclust()};
#'   fallback to \code{stats::hclust()}.
#'
#' @details
#' This function clusters **nodes** (not plots). The returned grouping can be fed
#' directly to \code{assign_releves(clusters = ...)}: we return
#' \code{attr(result, "clusters_for_assign")} as a **list of integer vectors**,
#' where each list element corresponds to one partition group \code{g1, g2, ...},
#' and contains the **node IDs** (integers) assigned to that group. If any nodes
#' were removed as duplicates (identical plot membership), we keep only the
#' rightmost occurrence.
#'
#' Group labels are always \code{"g1"}, \code{"g2"}, ... regardless of method.
#'
#' @return
#' A named character vector mapping node labels (e.g., \code{"c_123"}) to
#' partition labels (\code{"g1"}, \code{"g2"}, ...). The object has attributes:
#' \itemize{
#'   \item \code{"clusters_for_assign"}: \strong{list of integer vectors},
#'         one per group \code{g1..gK}, each containing the node IDs in that group.
#'   \item \code{"details"}: list with method, \code{k_partition}, kept node IDs,
#'         dedup info, and (for hclust) the \code{hclust} object and settings.
#' }
#'
#' @seealso \code{\link{assign_releves}}
#' @export
group_clusters <- function(
    x,
    phi_cut          = 0.30,
    k_partition      = 2L,
    method           = c("kmeans","hclust"),
    feature_scale    = TRUE,
    kmeans_seed      = NULL,
    dist_method      = "manhattan",
    clust_method     = "average",
    use_parallelDist = FALSE,
    use_fastcluster  = FALSE
) {
  ## --- helpers ---
  .stop_x <- function() {
    stop("x must include Cluster.height and Plot.cluster.", call. = FALSE)
  }
  .relabel_by_first_leaf <- function(hc, grp) {
    # Relabel groups by left-to-right position of their first leaf
    ord <- hc$order
    labs <- hc$labels
    pos <- match(labs, labs[ord])
    pos_first <- tapply(seq_along(labs), grp, function(ix) min(pos[ix]))
    old_ids <- as.integer(names(pos_first))
    new_ids <- seq_along(old_ids)[order(unlist(pos_first))]
    map <- stats::setNames(new_ids, old_ids)
    unname(map[as.character(grp)])
  }

  ## --- checks ---
  method <- match.arg(method)
  if (!is.list(x) || !all(c("Cluster.height","Plot.cluster") %in% names(x))) .stop_x()

  H  <- x$Cluster.height
  PC <- x$Plot.cluster
  if (is.null(colnames(PC))) stop("x$Plot.cluster must have node (cluster) column names like 'c_###'.")

  if (!is.numeric(phi_cut) || length(phi_cut) != 1L || !is.finite(phi_cut) || phi_cut < 0 || phi_cut > 1)
    stop("`phi_cut` must be a number in [0,1].", call. = FALSE)

  k_partition <- as.integer(k_partition)
  if (!is.finite(k_partition) || k_partition < 1L)
    stop("`k_partition` must be a positive integer.", call. = FALSE)

  ## --- 1) keep ALL nodes with H >= phi_cut (no topmost filtering) ---
  keep_nodes <- which(H >= phi_cut)
  if (!length(keep_nodes)) stop("No nodes meet phi_cut = ", phi_cut, call. = FALSE)
  kept_ids  <- keep_nodes
  kept_labs <- paste0("c_", kept_ids)

  # subset PC to kept nodes; build node × plot matrix
  if (!all(kept_labs %in% colnames(PC)))
    stop("Internal mismatch: some kept nodes not found in Plot.cluster columns.", call. = FALSE)

  M <- t(as.matrix(PC[, kept_labs, drop = FALSE]))  # rows = nodes, cols = plots
  storage.mode(M) <- "double"

  ## --- 3) deduplicate identical nodes (rows); keep the rightmost node ---
  # "Rightmost" = the node with the highest original column index in PC among duplicates.
  # We can achieve that by scanning in original kept_labs order and, for each key, keeping the LAST index.
  keys <- apply(M, 1L, function(r) paste0(as.integer(r), collapse = ""))  # binary signature per node-row
  # map key -> indices; keep last per key
  idx_by_key <- split(seq_along(keys), keys)
  keep_row_ix <- vapply(idx_by_key, function(ix) ix[length(ix)], integer(1))
  keep_row_ix <- sort(keep_row_ix)  # preserve overall order of last occurrences
  Mded <- M[keep_row_ix, , drop = FALSE]
  kept_ids_dedup  <- kept_ids[keep_row_ix]
  kept_labs_dedup <- kept_labs[keep_row_ix]

  if (nrow(Mded) < k_partition) {
    stop("After deduplication there are only ", nrow(Mded),
         " distinct nodes; cannot form k = ", k_partition, " groups.", call. = FALSE)
  }

  ## --- 4) feature scaling (optional) ---
  Muse <- if (isTRUE(feature_scale)) scale(Mded) else Mded

  ## --- 5) partition nodes ---
  if (method == "kmeans") {
    if (!is.null(kmeans_seed)) set.seed(as.integer(kmeans_seed))
    km <- stats::kmeans(Muse, centers = k_partition, iter.max = 100, nstart = 10)
    grp_ids <- km$cluster
    # label groups as g1..gK (ordered by numeric id of grp_ids)
    grp_labels <- paste0("g", grp_ids)

  } else {  # hclust
    if (use_parallelDist && !requireNamespace("parallelDist", quietly = TRUE)) {
      warning("parallelDist not available; falling back to stats::dist().")
      use_parallelDist <- FALSE
    }
    d <- if (use_parallelDist) parallelDist::parDist(Muse, method = dist_method)
    else stats::dist(Muse, method = dist_method)

    if (use_fastcluster && !requireNamespace("fastcluster", quietly = TRUE)) {
      warning("fastcluster not available; falling back to stats::hclust().")
      use_fastcluster <- FALSE
    }
    hc <- if (use_fastcluster) fastcluster::hclust(d, method = clust_method)
    else stats::hclust(d, method = clust_method)

    grp_raw <- stats::cutree(hc, k = k_partition)
    # Relabel by first-leaf order for readability, then convert to g1..gK
    grp_ord <- .relabel_by_first_leaf(hc, grp_raw)
    grp_labels <- paste0("g", grp_ord)
  }

  ## --- build outputs ---
  # Mapping node label -> group label (only for dedup-kept nodes)
  out <- stats::setNames(grp_labels, kept_labs_dedup)

  # clusters_for_assign: list of integer vectors (node IDs) per group g1..gK
  groups <- split(kept_ids_dedup, grp_labels)
  # order by g#, ensure stable order
  g_names <- paste0("g", sort(unique(as.integer(sub("^g", "", names(groups))))))
  clusters_for_assign <- lapply(g_names, function(gn) {
    if (gn %in% names(groups)) sort(groups[[gn]]) else integer(0)
  })
  names(clusters_for_assign) <- g_names

  # details
  details <- list(
    method        = method,
    k_partition   = k_partition,
    phi_cut       = phi_cut,
    feature_scale = feature_scale,
    dist_method   = if (method == "hclust") dist_method else NULL,
    clust_method  = if (method == "hclust") clust_method else NULL,
    kept_node_ids        = kept_ids,
    kept_node_ids_dedup  = kept_ids_dedup,
    kept_node_labels     = kept_labs,
    kept_node_labels_dedup = kept_labs_dedup
  )
  if (method == "hclust") details$hclust <- hc

  structure(out,
            clusters_for_assign = clusters_for_assign,
            details = details)
}
