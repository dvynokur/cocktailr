#' Select Cocktail clusters
#'
#' Selects Cocktail clusters from a \code{"cocktail"} object by first filtering
#' candidate clusters and then resolving nested clusters according to the
#' selected mode.
#'
#' Candidate clusters can be filtered by merge height (\code{min_phi}), species
#' group size (\code{min_k}), number of member plots (\code{min_n}), and score
#' (\code{min_score}). Optionally, selection can be restricted to clusters
#' contained within a broader cluster using \code{within_cluster}.
#'
#' The default score \code{"h_logk_logm"} is calculated as:
#' \deqn{
#' \mathrm{score}_c = h_c \cdot \log(k_c) \cdot \log(m_c),
#' }
#' where \code{h_c} is the Cocktail merge height, \code{k_c} is the number of
#' species in the cluster, and \code{m_c} is the Cocktail membership threshold.
#' This score is an operational ranking criterion, not a universal measure of
#' diagnostic quality. Use \code{score_method = "h_logk"} if the membership
#' threshold should not affect ranking but species-group size should still be
#' considered. Use \code{score_method = "h"} if selection should be based only
#' on merge height, i.e. the phi value at which each cluster was formed.
#'
#' Nesting is controlled by \code{mode}:
#' \itemize{
#'   \item \code{mode = "strict"}: no nesting allowed. After selecting a
#'   cluster, all its ancestors and descendants are excluded from further
#'   selection.
#'   \item \code{mode = "top"}: keep only the hierarchically highest eligible
#'   clusters by excluding candidates whose ancestors are also eligible.
#' }
#'
#' @param x A \code{"cocktail"} object, usually returned by
#'   \code{\link{cocktail_cluster}}, containing at least
#'   \code{Cluster.species}, \code{Cluster.height}, \code{Cluster.info}, and
#'   \code{Plot.cluster}.
#' @param clusters Optional vector of candidate clusters. Can be an integer
#'   vector of cluster IDs, e.g. \code{c(12, 27)}, or a character vector of
#'   labels, e.g. \code{c("c_12", "c_27")}. If \code{NULL}, all clusters are
#'   candidates.
#' @param within_cluster Optional broader cluster within which candidate
#'   clusters should be selected. Can be a single integer cluster ID or a
#'   character label, e.g. \code{"c_12"}. Candidate clusters are identified by
#'   topological species-set containment in \code{x$Cluster.species}. The
#'   \code{within_cluster} itself is not returned as a candidate.
#' @param min_phi Numeric. Minimum merge height, i.e. the phi value at which the
#'   Cocktail cluster was formed. Default \code{0.2}.
#' @param min_k Integer. Minimum number of species in a cluster. Default
#'   \code{1L}.
#' @param min_n Integer. Minimum number of plots satisfying the
#'   cluster-specific Cocktail membership rule. Default \code{1L}.
#' @param min_score Numeric. Minimum selection score required for a cluster to
#'   be eligible. Default \code{0}.
#' @param score_method Character. Method used to rank candidate clusters:
#'   \itemize{
#'     \item \code{"h_logk_logm"}: \code{h * log(k) * log(m)};
#'     \item \code{"h_logk"}: \code{h * log(k)};
#'     \item \code{"h"}: merge height only.
#'   }
#' @param mode Character. How nested clusters are handled during selection.
#'   Either \code{"strict"} or \code{"top"}.
#' @param return Character. Type of output:
#'   \itemize{
#'     \item \code{"labels"}: character labels like \code{"c_12"};
#'     \item \code{"ids"}: integer cluster IDs;
#'     \item \code{"table"}: a data frame with selected clusters and their
#'     \code{id}, \code{h}, \code{k}, \code{m}, \code{n}, and \code{score}.
#'   }
#'
#' @return Depending on \code{return}, a character vector, an integer vector, or
#'   a data frame with selected Cocktail clusters.
#'
#' @export

select_clusters <- function(
    x,
    clusters = NULL,
    within_cluster = NULL,
    min_phi = 0.2,
    min_k = 1L,
    min_n = 1L,
    min_score = 0,
    score_method = c("h_logk_logm", "h_logk", "h"),
    mode = c("strict", "top"),
    return = c("labels", "ids", "table")
) {
  return <- match.arg(return)
  mode <- match.arg(mode)
  score_method <- match.arg(score_method)

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
  if (!"Plot.cluster" %in% names(x) || is.null(x$Plot.cluster)) {
    stop("`x$Plot.cluster` is missing; cannot compute the number of member plots `n`.")
  }

  CS <- x$Cluster.species
  H <- x$Cluster.height
  KI <- x$Cluster.info
  PC <- x$Plot.cluster

  if (!is.matrix(CS)) stop("`x$Cluster.species` must be a matrix.")
  n_nodes <- nrow(CS)

  if (length(H) < n_nodes) {
    stop("`x$Cluster.height` must have length >= nrow(x$Cluster.species).")
  }
  if (ncol(PC) < n_nodes) {
    stop("`x$Plot.cluster` must have at least nrow(x$Cluster.species) columns.")
  }

  if (!is.numeric(min_phi) || length(min_phi) != 1L || is.na(min_phi)) {
    stop("`min_phi` must be a single numeric value.")
  }
  if (!is.numeric(min_k) || length(min_k) != 1L || is.na(min_k)) {
    stop("`min_k` must be a single integer-like value.")
  }
  if (!is.numeric(min_n) || length(min_n) != 1L || is.na(min_n)) {
    stop("`min_n` must be a single integer-like value.")
  }
  if (!is.numeric(min_score) || length(min_score) != 1L || is.na(min_score)) {
    stop("`min_score` must be a single numeric value.")
  }

  min_k <- as.integer(min_k)
  min_n <- as.integer(min_n)

  ## ---- helper functions ---------------------------------------------------
  parse_cluster_ids <- function(z, arg_name = "clusters") {
    if (is.null(z)) return(integer(0))
    if (is.list(z)) z <- unlist(z, use.names = FALSE)

    if (is.character(z)) {
      out <- suppressWarnings(as.integer(sub("^c_", "", z)))
    } else {
      out <- suppressWarnings(as.integer(z))
    }

    out <- out[is.finite(out) & out > 0L & out <= n_nodes]
    out <- sort(unique(out))

    if (!length(out)) {
      stop("No valid cluster IDs found in `", arg_name, "` after filtering to 1..", n_nodes, ".")
    }
    out
  }

  species_set <- function(node) which(CS[node, ] > 0L)

  set_is_subset <- function(a, b) {
    Sa <- species_set(a)
    Sb <- species_set(b)
    if (length(Sa) == 0L) return(TRUE)
    if (length(Sa) > length(Sb)) return(FALSE)
    all(Sa %in% Sb)
  }

  set_is_strict_subset <- function(a, b) {
    set_is_subset(a, b) && !set_is_subset(b, a)
  }

  compute_score <- function(h, k, m, method) {
    score <- switch(
      method,
      h_logk_logm = h * log(k) * log(m),
      h_logk      = h * log(k),
      h           = h
    )
    score[!is.finite(score)] <- 0
    score
  }

  ## ---- parse clusters argument into candidate IDs -------------------------
  if (missing(clusters) || is.null(clusters)) {
    ids <- seq_len(n_nodes)
  } else {
    ids <- parse_cluster_ids(clusters, "clusters")
  }

  ## ---- optionally restrict to descendants of within_cluster ---------------
  if (!is.null(within_cluster)) {
    parent <- parse_cluster_ids(within_cluster, "within_cluster")
    if (length(parent) != 1L) {
      stop("`within_cluster` must identify exactly one cluster.")
    }
    parent <- parent[[1L]]

    ids <- ids[vapply(ids, function(i) set_is_strict_subset(i, parent), logical(1))]
    if (!length(ids)) {
      stop("No candidate clusters are strict descendants of `within_cluster = c_", parent, "`.")
    }
  }

  ## ---- compute candidate statistics ---------------------------------------
  h <- as.numeric(H[ids])
  k <- as.numeric(KI[ids, "k"])
  m <- as.numeric(KI[ids, "m"])
  n <- as.numeric(Matrix::colSums(PC[, ids, drop = FALSE] > 0))

  h[!is.finite(h)] <- 0
  k[!is.finite(k) | k < 0] <- 0
  m[!is.finite(m) | m < 0] <- 0
  n[!is.finite(n) | n < 0] <- 0

  score <- compute_score(h, k, m, score_method)

  cand_tbl <- data.frame(
    cluster = paste0("c_", ids),
    id = ids,
    h = h,
    k = k,
    m = m,
    n = n,
    score = score,
    stringsAsFactors = FALSE
  )

  ## ---- filtering ----------------------------------------------------------
  keep <- cand_tbl$h >= min_phi &
    cand_tbl$k >= min_k &
    cand_tbl$n >= min_n &
    cand_tbl$score >= min_score

  cand_tbl <- cand_tbl[keep, , drop = FALSE]

  if (!nrow(cand_tbl)) {
    stop(
      "No clusters remain after applying selection filters: ",
      "min_phi = ", min_phi, ", min_k = ", min_k,
      ", min_n = ", min_n, ", min_score = ", min_score, "."
    )
  }

  ## ---- ordering: score -> h -> k -> n -> m -> smaller id ------------------
  ord <- order(
    -cand_tbl$score,
    -cand_tbl$h,
    -cand_tbl$k,
    -cand_tbl$n,
    -cand_tbl$m,
    cand_tbl$id,
    na.last = NA
  )
  cand_tbl <- cand_tbl[ord, , drop = FALSE]
  cand <- cand_tbl$id

  ## ---- selection ----------------------------------------------------------
  if (mode == "strict") {
    selected <- integer(0)

    for (node in cand) {
      if (!length(selected)) {
        selected <- c(selected, node)
        next
      }

      node_is_descendant <- any(vapply(selected, function(s) set_is_subset(node, s), logical(1)))
      node_is_ancestor <- any(vapply(selected, function(s) set_is_subset(s, node), logical(1)))

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
        if (set_is_strict_subset(selected[j], selected[i])) {
          drop[j] <- TRUE
        }
      }
    }

    selected <- selected[!drop]
  }

  ## ---- return -------------------------------------------------------------
  if (return == "ids") {
    return(selected)
  }

  if (return == "labels") {
    return(paste0("c_", selected))
  }

  out <- cand_tbl[match(selected, cand_tbl$id), , drop = FALSE]
  rownames(out) <- NULL
  out
}
