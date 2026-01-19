#' Select Cocktail clusters using a combined score based on merge phi, size, and m-threshold
#'
#' @description
#' Selects internal Cocktail clusters (nodes) by ranking them with a combined score
#' based on the Cocktail merge height \code{h} (phi), the number of species \code{k},
#' and the Cocktail membership threshold \code{m}.
#'
#' The scoring function is:
#' \deqn{
#' \mathrm{score}(k) = h_k \cdot \log(k_k) \cdot \log(m_k),
#' }
#' where:
#' \itemize{
#'   \item \code{h_k} is \code{x$Cluster.height[k]} (merge phi),
#'   \item \code{k_k} is \code{x$Cluster.info[k,"k"]} (cluster size),
#'   \item \code{m_k} is \code{x$Cluster.info[k,"m"]} (Cocktail membership threshold).
#' }
#'
#' Clusters are selected strongest-first, using a greedy procedure in which eligible
#' clusters are considered in descending score order (ties broken by larger \code{k},
#' then smaller node ID).
#'
#' Nesting is controlled by \code{mode}:
#' \itemize{
#'   \item \code{mode = "strict"} (default): \strong{no nesting allowed}.
#'     After selecting a cluster, \emph{all} of its strict ancestors and strict
#'     descendants are excluded. Therefore, the final selection contains no nested
#'     clusters (no strict subset/superset relationships between selected clusters
#'     based on topological species sets from \code{x$Cluster.species}).
#'   \item \code{mode = "top"}: \strong{hierarchically highest clusters among eligible}.
#'     All eligible clusters (\code{score >= min_score}) are first collected.
#'     The final selection then keeps only those clusters that are not strict
#'     descendants of another eligible cluster. In other words, if a cluster has an
#'     eligible ancestor, it is removed and the ancestor is kept.
#' }
#'
#' @param x A \code{"cocktail"} object (result of \code{\link{cocktail_cluster}}),
#'   containing at least \code{Cluster.species}, \code{Cluster.height}, and
#'   \code{Cluster.info} with columns \code{"k"} and \code{"m"}.
#'
#' @param clusters Optional cluster identifiers (nodes) to consider. Can be a numeric
#'   vector of node indices (e.g. \code{c(12, 27)}) or a character vector of labels
#'   (e.g. \code{c("c_12", "c_27")}). Each element refers to a single internal node.
#'   If \code{NULL} or missing, all internal nodes \code{1:nrow(x$Cluster.species)}
#'   are candidates.
#'
#' @param min_phi Numeric scalar; minimum merge phi required to keep a node.
#'   Default \code{0.2}. Nodes with \code{Cluster.height < min_phi} are dropped.
#'
#' @param min_k Integer; minimum number of species \code{k} required for a cluster
#'   to be eligible. Default \code{1}. This can be used to filter out very small
#'   clusters.
#'
#' @param min_score Numeric scalar; minimum score required for a cluster to be eligible.
#'   Default \code{1}. Set to \code{0} to allow the full selection (no score filtering).
#'
#' @param mode Character; nesting rule used when constructing the final selection.
#'   \itemize{
#'     \item \code{"strict"} (default): exclude both ancestors and descendants of
#'       already selected clusters (no nested clusters are returned).
#'     \item \code{"top"}: keep only hierarchically highest clusters among the
#'       eligible set (drop descendants if an eligible ancestor exists).
#'   }
#'
#' @param return What to return:
#'   \itemize{
#'     \item \code{"labels"} (default): character labels like \code{"c_12"}.
#'     \item \code{"ids"}: integer node IDs.
#'     \item \code{"table"}: a data frame with selected clusters and their \code{h,k,m,score}.
#'   }
#'
#' @return Depending on \code{return}:
#' \itemize{
#'   \item \code{"labels"}: character vector of selected cluster labels.
#'   \item \code{"ids"}: integer vector of selected node IDs.
#'   \item \code{"table"}: data frame of selected clusters sorted by decreasing score.
#' }
#'
#' @export
select_clusters <- function(
    x,
    clusters = NULL,
    min_phi = 0.2,
    min_k = 1L,
    min_score = 1,
    mode = c("strict", "top"),
    return = c("labels", "ids", "table")
) {
  return <- match.arg(return)
  mode   <- match.arg(mode)

  ## ---- basic checks -------------------------------------------------------
  if (!is.list(x) || !"Cluster.species" %in% names(x)) {
    stop("`x` must be a Cocktail object with a `Cluster.species` component.")
  }
  if (!"Cluster.height" %in% names(x) || is.null(x$Cluster.height)) {
    stop("`x$Cluster.height` is missing; cannot compute h-based score.")
  }
  if (!"Cluster.info" %in% names(x) || is.null(x$Cluster.info)) {
    stop("`x$Cluster.info` is missing; cannot compute k and m.")
  }
  if (!all(c("k", "m") %in% colnames(x$Cluster.info))) {
    stop("`x$Cluster.info` must contain columns 'k' and 'm'.")
  }

  CS <- x$Cluster.species
  H  <- x$Cluster.height
  KI <- x$Cluster.info

  if (!is.matrix(CS)) stop("`x$Cluster.species` must be a matrix.")
  n_nodes <- nrow(CS)

  if (!is.numeric(min_phi) || length(min_phi) != 1L || is.na(min_phi)) {
    stop("`min_phi` must be a single numeric value.")
  }

  if (!is.numeric(min_k) || length(min_k) != 1L || is.na(min_k)) {
    stop("`min_k` must be a single integer-like value.")
  }
  min_k <- as.integer(min_k)

  if (!is.numeric(min_score) || length(min_score) != 1L || is.na(min_score)) {
    stop("`min_score` must be a single numeric value.")
  }

  if (length(H) < n_nodes) {
    stop("`x$Cluster.height` must have length >= nrow(x$Cluster.species).")
  }

  ## ---- parse clusters argument into candidate IDs -------------------------
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

  ## ---- filter by min_phi --------------------------------------------------
  ids <- ids[H[ids] >= min_phi]
  if (!length(ids)) {
    stop("No clusters remain after applying `min_phi = ", min_phi, "`.")
  }

  ## ---- filter by min_k ----------------------------------------------------
  k_all <- as.numeric(KI[ids, "k"])
  keep_k <- is.finite(k_all) & (k_all >= min_k)
  ids <- ids[keep_k]
  if (!length(ids)) {
    stop("No clusters remain after applying `min_k = ", min_k, "`.")
  }

  ## ---- compute score = h * log(k) * log(m) -------------------------------
  h <- as.numeric(H[ids])
  k <- as.numeric(KI[ids, "k"])
  m <- as.numeric(KI[ids, "m"])

  # guard against weird values
  h[!is.finite(h)] <- 0
  k[!is.finite(k) | k < 0] <- 0
  m[!is.finite(m) | m < 0] <- 0

  score <- h * log(k) * log(m)

  ## ---- filter by min_score ------------------------------------------------
  keep_score <- is.finite(score) & (score >= min_score)
  ids   <- ids[keep_score]
  h     <- h[keep_score]
  k     <- k[keep_score]
  m     <- m[keep_score]
  score <- score[keep_score]

  if (!length(ids)) {
    stop("No clusters remain after applying `min_score = ", min_score, "`.")
  }

  ## ---- tie-breaking: prefer larger k, then lower node id -------------------
  ord <- order(score, k, ids, decreasing = c(TRUE, TRUE, FALSE), na.last = NA)
  cand <- ids[ord]

  ## ---- precompute species sets for nesting checks --------------------------
  sp_sets <- lapply(cand, function(node) which(CS[node, ] > 0L))
  names(sp_sets) <- as.character(cand)

  is_subset <- function(a, b) {
    Sa <- sp_sets[[as.character(a)]]
    Sb <- sp_sets[[as.character(b)]]
    if (length(Sa) == 0L) return(TRUE)
    if (length(Sa) > length(Sb)) return(FALSE)
    all(Sa %in% Sb)
  }

  ## ---- selection ----------------------------------------------------------
  if (mode == "strict") {
    # greedy selection; exclude ancestors and descendants of already selected nodes
    selected <- integer(0)

    for (node in cand) {
      if (!length(selected)) {
        selected <- c(selected, node)
        next
      }

      node_is_descendant <- any(vapply(selected, function(s) is_subset(node, s), logical(1)))
      node_is_ancestor   <- any(vapply(selected, function(s) is_subset(s, node), logical(1)))

      if (node_is_descendant || node_is_ancestor) next
      selected <- c(selected, node)
    }

  } else {
    # mode == "top": keep only hierarchically highest clusters among eligible candidates
    selected <- cand
    drop <- logical(length(selected))

    for (i in seq_along(selected)) {
      if (drop[i]) next
      for (j in seq_along(selected)) {
        if (i == j || drop[j]) next
        # if selected[j] is a strict descendant of selected[i], drop it
        if (is_subset(selected[j], selected[i]) && !is_subset(selected[i], selected[j])) {
          drop[j] <- TRUE
        }
      }
    }

    selected <- selected[!drop]
  }

  ## ---- return --------------------------------------------------------------
  if (return == "ids") {
    return(selected)
  }

  if (return == "labels") {
    return(paste0("c_", selected))
  }

  # return == "table"
  out <- data.frame(
    cluster = paste0("c_", ids),
    h       = h,
    k       = k,
    m       = m,
    score   = score,
    stringsAsFactors = FALSE
  )
  out <- out[match(paste0("c_", selected), out$cluster), , drop = FALSE]
  out
}
