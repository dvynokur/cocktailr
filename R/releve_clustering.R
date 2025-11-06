#' Hierarchical clustering of relevés based on shared Cocktail clusters
#'
#' @description
#' Builds a hierarchical clustering (`hclust`) of plots from shared cluster
#' memberships (`Plot.cluster`). Steps:
#'  1) Keep clusters with φ >= `phi_cut`.
#'  2) Remove redundant cluster columns with identical plot membership (keep last).
#'  3) Compute distances (optionally with parallelDist).
#'  4) Agglomerate (optionally with fastcluster).
#'
#' @param x A Cocktail object (from `cocktail_cluster_mod()` or `cocktail_cluster()`).
#' @param phi_cut Numeric. Keep only clusters with height ≥ `phi_cut`. Default = 0.2.
#' @param dist_method Character. Distance for plots (default "manhattan").
#' @param clust_method Character. Linkage for hierarchical clustering (default "average").
#' @param dist_engine  "base" (stats::dist) or "parallelDist" (parallelDist::parDist). Default "base".
#' @param hclust_engine "stats" (stats::hclust) or "fastcluster" (fastcluster::hclust). Default "stats".
#' @param par_cores Integer or NULL. Cores for parallelDist; NULL lets parallelDist choose.
#'
#' @return A list of class `"releve_hclust"` with:
#'   - `hclust`         : the hclust object
#'   - `dist`           : the distance object used
#'   - `kept_clusters`  : indices of clusters retained (after φ cut & dedup)
#'   - `plot_labels`    : plot names
#'   - `engines`        : list(dist_engine, hclust_engine)
#' @export
releve_clustering <- function(
    x,
    phi_cut = 0.2,
    dist_method = "manhattan",
    clust_method = "average",
    dist_engine = c("base", "parallelDist"),
    hclust_engine = c("stats", "fastcluster"),
    par_cores = NULL
) {
  # --- validate input ---
  need <- c("Cluster.species", "Plot.cluster", "Cluster.height")
  if (!is.list(x) || !all(need %in% names(x))) {
    stop("x must be a valid Cocktail object (from cocktail_cluster_mod or cocktail_cluster).")
  }
  dist_engine  <- match.arg(dist_engine)
  hclust_engine <- match.arg(hclust_engine)

  PC <- x$Plot.cluster
  H  <- x$Cluster.height
  plot_names <- rownames(PC); if (is.null(plot_names)) plot_names <- as.character(seq_len(nrow(PC)))

  # --- 1) φ cut ---
  keep <- which(H >= phi_cut)
  if (length(keep) == 0L)
    stop("No clusters remain after applying phi_cut = ", phi_cut)

  PCk <- PC[, keep, drop = FALSE]

  # --- 2) remove redundant cluster columns by identical plot membership (keep last) ---
  keys <- apply(PCk, 2L, paste0, collapse = "")
  last_keep <- unlist(lapply(split(seq_along(keys), keys), tail, 1))
  ord_last  <- sort(last_keep)
  keep <- keep[ord_last]
  PCk  <- PCk[, ord_last, drop = FALSE]

  if (ncol(PCk) == 0L)
    stop("No informative clusters remain after filtering (all redundant).")

  # --- 3) distance matrix (optionally parallel) ---
  # graceful fallback if package not available
  if (dist_engine == "parallelDist") {
    if (!requireNamespace("parallelDist", quietly = TRUE)) {
      warning("parallelDist not installed; falling back to stats::dist().")
      dist_engine <- "base"
    }
  }

  d <- switch(
    dist_engine,
    base = stats::dist(PCk, method = dist_method),
    parallelDist = {
      # parallelDist::parDist returns an object of class 'dist'
      if (is.null(par_cores)) {
        parallelDist::parDist(PCk, method = dist_method)
      } else {
        parallelDist::parDist(PCk, method = dist_method, threads = par_cores)
      }
    }
  )
  # attach labels (some downstream code may rely on them)
  attr(d, "Labels") <- plot_names

  # --- 4) hierarchical clustering (optionally fastcluster) ---
  if (hclust_engine == "fastcluster") {
    if (!requireNamespace("fastcluster", quietly = TRUE)) {
      warning("fastcluster not installed; falling back to stats::hclust().")
      hclust_engine <- "stats"
    }
  }

  hc <- switch(
    hclust_engine,
    stats = stats::hclust(d, method = clust_method),
    fastcluster = fastcluster::hclust(d, method = clust_method)
  )

  structure(
    list(
      hclust = hc,
      dist = d,
      kept_clusters = keep,
      plot_labels = plot_names,
      engines = list(dist = dist_engine, hclust = hclust_engine)
    ),
    class = "releve_hclust"
  )
}
