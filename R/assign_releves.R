#' Assign relevés (plots) to Cocktail groups or partition them
#'
#' @description
#' Assign plots using the Cocktail tree (`x`) in several ways:
#'
#' - **"topology_only"** *(default)* — Use topology from `x$Plot.cluster`
#'   at a given `phi_cut` (or a user-specified set of parent clusters).
#'   Each plot is assigned to the **unique** group (if any) where it is a member
#'   at that cut; ties or no-membership → `"not assigned"`.
#'
#' - **"topology_cover"** — Like *topology_only*, but when `vegmatrix` contains cover
#'   values, score parent groups by **sum of covers** of their member species in each plot,
#'   restricted to plots that meet the group’s m-rule in `x$Plot.cluster`.
#'   Scores are **per-plot normalised** (divide by the row sum across candidate groups)
#'   before taking a unique maximum. Plots failing `min_cover` for a group are ineligible.
#'
#' - **"cover_strict"** — Ignore topology membership and score each plot–group pair as the
#'   **sum of transformed covers** across the group’s species. Unique maximum wins; otherwise
#'   `"not assigned"`. Uses only parent groups at `phi_cut` (or `clusters` if given).
#'
#' - **"cover_fuzzy"** — Compute a species × node \eqn{\phi} matrix from `x$Plot.cluster`,
#'   restrict to parent groups at `phi_cut` (or `clusters`), set non-positive \eqn{\phi} to 0,
#'   then score plots by (transformed cover) × \eqn{\phi} weights. Unique maximum wins.
#'
#' - **"partition_kmeans"** — Build a feature matrix from `x$Plot.cluster` filtered to
#'   parent groups at `phi_cut` (or `clusters`), **deduplicate identical columns**, optionally
#'   **standardize features** (`feature_scale = TRUE`), then do k-means into `k_partition` groups.
#'   For large vegetation tables, this can be slow.
#'
#' - **"partition_hclust"** — Same features as above; compute distances (default
#'   Manhattan), then hierarchical clustering (`clust_method`, default UPGMA/"average")
#'   and cut into `k_partition` groups. **Clusters are relabeled** in left-to-right
#'   dendrogram order for readability. For large vegetation tables, this can be slow.
#'
#' @param x A Cocktail result (list) with at least `Cluster.species`, `Cluster.merged`,
#'   `Cluster.height`, and `Plot.cluster` (from `cocktail_cluster()`).
#' @param vegmatrix Optional numeric matrix/data.frame of **covers** with plots in rows
#'   and species in columns (needed for *_cover* strategies). Missing values are treated as 0.
#' @param strategy One of
#'   `c("topology_only","topology_cover","cover_strict","cover_fuzzy",
#'      "partition_kmeans","partition_hclust")`.
#' @param phi_cut Numeric between 0 and 1 (inclusive). Cut level used to choose
#'   **parent clusters** (unless `clusters` is supplied).
#' @param clusters Optional vector of cluster IDs to use instead of `phi_cut`.
#'   Accepts character labels like `"c_1333"` or integer IDs like `1333`.
#' @param cover_transform Cover transform used by cover-based strategies:
#'   `c("sqrt","none","log1p")`. Default `"sqrt"`.
#' @param min_cover Numeric percent between 0 and 100 (inclusive). A plot is eligible
#'   for a group only if the **raw (pre-transform)** sum of that group’s species covers
#'   in the plot is at least this value. Default 0.
#' @param k_partition Positive integer, required for partitioning strategies.
#' @param kmeans_seed Optional scalar to set the RNG seed for k-means reproducibility.
#' @param dist_method Distance for `"partition_hclust"` (default `"manhattan"`).
#' @param clust_method Linkage for `"partition_hclust"` (default `"average"`).
#' @param use_parallelDist Logical. If `TRUE`, try `parallelDist::parDist()` for distance.
#'   Falls back to `stats::dist()` if the package is unavailable.
#' @param use_fastcluster Logical. If `TRUE`, try `fastcluster::hclust()` for agglomeration.
#'   Falls back to `stats::hclust()` if the package is unavailable.
#' @param feature_scale Logical; if `TRUE` (default), standardize features
#'   (columns) for partitioning strategies to equalize scale.
#' @param min_group_size Integer at least 1. After assignment/partitioning, any group with fewer
#'   than this many plots is relabeled `"not assigned"`. Default 1 (no collapsing).
#'
#' @details
#' **Parent clusters** are the “topmost” nodes with height at least `phi_cut`
#' (no ancestor also at or above the cut). If `clusters` is supplied, it must
#' identify valid Cocktail node IDs/labels; mixed types are allowed (e.g., `c("c_12", 25)`).
#'
#' For cover-based scoring, species columns of `vegmatrix` are aligned to
#' `colnames(x$Cluster.species)`. Non-overlapping species are ignored. If your covers are
#' on **ordinal** scales (e.g., 1–5, Braun–Blanquet codes), results are computed but a
#' message is emitted recommending percentage covers for best performance.
#'
#' **Runtime note:** partitioning on **large vegetation tables** (many plots and/or many
#' clusters kept at the cut) can take a long time due to distance and clustering steps.
#'
#' @return
#' A named character vector of group labels per plot (or `"not assigned"`). The return
#' value has an attribute `"details"` with:
#' \itemize{
#'   \item `strategy`, `phi_cut`, `clusters_used`, `feature_scale`
#'   \item `scores` (for scoring strategies) or `partition` (data.frame with cluster IDs)
#'   \item `parents` (parent node IDs used)
#'   \item `min_cover`, `cover_transform`
#'   \item `min_group_size`, `collapsed_groups`
#' }
#'
#' @seealso \code{\link{clusters_at_cut}}
#' @export

assign_releves <- function(
    x,
    vegmatrix        = NULL,
    strategy         = c("topology_only","topology_cover","cover_strict","cover_fuzzy",
                         "partition_kmeans","partition_hclust"),
    phi_cut          = 0.30,
    clusters         = NULL,
    cover_transform  = c("sqrt","none","log1p"),
    min_cover        = 0,
    k_partition      = NULL,
    kmeans_seed      = NULL,
    dist_method      = "manhattan",
    clust_method     = "average",
    use_parallelDist = FALSE,
    use_fastcluster  = FALSE,
    feature_scale    = TRUE,
    min_group_size   = 1
) {
  ## ---- helpers ----
  .is_cocktail <- function(obj) {
    is.list(obj) && all(c("Cluster.species","Cluster.merged","Cluster.height","Plot.cluster") %in% names(obj))
  }
  .cover_transform <- function(xx, how) {
    if (how == "sqrt") sqrt(pmax(xx, 0)) else if (how == "log1p") log1p(pmax(xx, 0)) else xx
  }
  .parents_at_cut <- function(CM, H, cut) {
    idx <- which(H >= cut); if (!length(idx)) return(integer(0))
    kids <- CM[idx, , drop = FALSE]
    children <- unique(as.integer(kids[kids > 0]))
    sort(setdiff(idx, intersect(idx, children)))
  }
  .parse_clusters_arg <- function(v) {
    if (is.null(v) || !length(v)) return(NULL)
    lab_like <- suppressWarnings(grepl("^c_\\d+$", v))
    ids <- rep(NA_integer_, length(v))
    ids[lab_like] <- as.integer(sub("^c_", "", v[lab_like]))
    ids[!lab_like] <- suppressWarnings(as.integer(v[!lab_like]))
    ids
  }
  .dedup_cols_by_key <- function(M) {
    if (!ncol(M)) return(list(M = M, keep = integer(0)))
    keys <- apply(M, 2L, function(col) paste0(as.integer(round(as.numeric(col), 12)), collapse = ""))
    last_keep <- vapply(split(seq_along(keys), keys), function(ix) ix[length(ix)], integer(1))
    ord <- sort(last_keep)
    list(M = M[, ord, drop = FALSE], keep = ord)
  }
  .relabel_by_first_leaf <- function(hc, grp) {
    ord <- hc$order
    labs <- hc$labels
    # map labels to integer positions
    pos <- match(labs, labs[ord])
    pos_first <- tapply(seq_along(labs), grp, function(ix) min(pos[ix]))
    old_ids <- as.integer(names(pos_first))
    new_ids <- seq_along(old_ids)[order(unlist(pos_first))]
    map <- stats::setNames(new_ids, old_ids)
    unname(map[as.character(grp)])
  }

  ## ---- args / checks ----
  strategy        <- match.arg(strategy)
  cover_transform <- match.arg(cover_transform)

  if (!.is_cocktail(x)) {
    stop("x must be a Cocktail object with Cluster.species, Cluster.merged, Cluster.height, and Plot.cluster.")
  }
  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  PC <- x$Plot.cluster
  spp_all <- colnames(CS)
  plot_names <- rownames(PC); if (is.null(plot_names)) plot_names <- paste0("plot_", seq_len(nrow(PC)))

  # choose parent nodes
  parent_ids <- if (!is.null(clusters)) {
    ids <- .parse_clusters_arg(clusters)
    if (any(is.na(ids))) stop("`clusters` must be labels like 'c_123' or integer node IDs.")
    ids
  } else {
    .parents_at_cut(CM, H, phi_cut)
  }
  if (!length(parent_ids)) stop("No parent clusters available (check phi_cut or clusters).")
  parent_labs <- paste0("c_", parent_ids)

  # covers (needed for *_cover)
  need_cover <- strategy %in% c("topology_cover","cover_strict","cover_fuzzy")
  if (need_cover) {
    if (is.null(vegmatrix)) stop("`vegmatrix` (covers) is required for strategy = ", strategy, ".")
    vm <- as.matrix(vegmatrix); vm[is.na(vm)] <- 0
    storage.mode(vm) <- "double"
    if (is.null(colnames(vm))) stop("`vegmatrix` must have species column names.")
    if (is.null(rownames(vm))) rownames(vm) <- plot_names
    # align
    common <- intersect(colnames(vm), spp_all)
    if (!length(common)) stop("No overlapping species between covers and Cocktail species.")
    vm <- vm[, common, drop = FALSE]
    # quick diagnostics
    if (all(vm %in% c(0,1))) {
      message("Covers appear binary; cover-based strategies will behave similarly to presence.")
    } else if (length(sort(unique(as.numeric(vm)))) <= 10) {
      message("Covers look ordinal/limited in range; percentage covers are recommended for best results.")
    }
  } else {
    vm <- NULL
  }

  # hygiene
  min_cover <- max(0, min(100, as.numeric(min_cover)))
  min_group_size <- as.integer(min_group_size)
  if (!is.finite(min_group_size) || min_group_size < 1) min_group_size <- 1L

  ## ---- group species sets and membership ----
  groups_species <- lapply(parent_ids, function(i) spp_all[CS[i, ] == 1L])
  names(groups_species) <- parent_labs

  keep_cols <- intersect(colnames(PC), parent_labs)
  if (!length(keep_cols)) stop("Internal mismatch: no Plot.cluster columns for chosen parents.")
  Gcut <- as.matrix(PC[, keep_cols, drop = FALSE]); storage.mode(Gcut) <- "double"

  assigned <- rep("not assigned", nrow(Gcut)); names(assigned) <- plot_names
  details <- list(strategy = strategy, phi_cut = phi_cut, clusters_used = parent_labs,
                  min_cover = min_cover, cover_transform = cover_transform,
                  feature_scale = feature_scale,
                  min_group_size = min_group_size)

  ## ---- strategies ----
  if (strategy == "topology_only") {
    hits <- Gcut > 0
    row_hits <- rowSums(hits)
    winner <- max.col(hits, ties.method = "first")
    ok <- row_hits == 1
    assigned[ok] <- colnames(Gcut)[winner[ok]]
    details$scores <- NULL

  } else if (strategy == "topology_cover") {
    S_raw <- matrix(0, nrow = nrow(Gcut), ncol = ncol(Gcut),
                    dimnames = list(rownames(vm), colnames(Gcut)))
    for (g in seq_along(groups_species)) {
      sp <- intersect(groups_species[[g]], colnames(vm))
      if (length(sp)) S_raw[, g] <- rowSums(vm[, sp, drop = FALSE])
    }
    eligible <- (Gcut > 0) & (S_raw >= min_cover)
    S <- S_raw
    S[!eligible] <- -Inf
    row_tot <- rowSums(pmax(S, 0), na.rm = TRUE)
    scale_vec <- ifelse(row_tot > 0, 1 / row_tot, 0)
    S_scaled <- S * scale_vec
    winner <- max.col(S_scaled, ties.method = "first")
    best <- S_scaled[cbind(seq_len(nrow(S_scaled)), winner)]
    ties <- apply(S_scaled, 1, function(v) sum(v == max(v)))
    ok <- is.finite(best) & (ties == 1)
    assigned[ok] <- colnames(S_scaled)[winner[ok]]
    details$scores <- S_scaled

  } else if (strategy == "cover_strict") {
    W <- .cover_transform(vm, cover_transform)
    S_raw <- matrix(0, nrow = nrow(W), ncol = length(groups_species),
                    dimnames = list(rownames(W), names(groups_species)))
    for (g in seq_along(groups_species)) {
      sp <- intersect(groups_species[[g]], colnames(W))
      if (length(sp)) S_raw[, g] <- rowSums(W[, sp, drop = FALSE])
    }
    S_elig <- matrix(0, nrow = nrow(vm), ncol = length(groups_species),
                     dimnames = list(rownames(vm), names(groups_species)))
    for (g in seq_along(groups_species)) {
      sp <- intersect(groups_species[[g]], colnames(vm))
      if (length(sp)) S_elig[, g] <- rowSums(vm[, sp, drop = FALSE])
    }
    S <- S_raw
    S[S_elig < min_cover] <- -Inf
    winner <- max.col(S, ties.method = "first")
    best <- S[cbind(seq_len(nrow(S)), winner)]
    ties <- apply(S, 1, function(v) sum(v == max(v)))
    ok <- is.finite(best) & (ties == 1)
    assigned[ok] <- colnames(S)[winner[ok]]
    details$scores <- S

  } else if (strategy == "cover_fuzzy") {
    X <- (vm > 0); storage.mode(X) <- "double"
    N <- nrow(X)
    G <- as.matrix(PC > 0); storage.mode(G) <- "double"
    if (nrow(G) != N) stop("Internal mismatch: Plot.cluster has different number of plots.")
    a  <- crossprod(X, G)
    p  <- colSums(X)
    g1 <- colSums(G)
    b <- matrix(p,  nrow = ncol(X), ncol = ncol(G)) - a
    c <- matrix(g1, nrow = ncol(X), ncol = ncol(G), byrow = TRUE) - a
    d <- (N - matrix(g1, nrow = ncol(X), ncol = ncol(G), byrow = TRUE)) - b
    den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
    Phi_all <- (a * d - b * c) / den
    Phi_all[!is.finite(den) | den <= 0] <- 0
    colnames(Phi_all) <- colnames(PC)
    rownames(Phi_all) <- colnames(X)
    Phi <- Phi_all[, parent_labs, drop = FALSE]
    Phi[Phi <= 0] <- 0
    W <- .cover_transform(vm, cover_transform)
    S <- as.matrix(W %*% Phi)
    S_elig <- matrix(0, nrow = nrow(vm), ncol = length(groups_species),
                     dimnames = list(rownames(vm), names(groups_species)))
    for (g in seq_along(groups_species)) {
      sp <- intersect(groups_species[[g]], colnames(vm))
      if (length(sp)) S_elig[, g] <- rowSums(vm[, sp, drop = FALSE])
    }
    S[S_elig < min_cover] <- -Inf
    winner <- max.col(S, ties.method = "first")
    best <- S[cbind(seq_len(nrow(S)), winner)]
    ties <- apply(S, 1, function(v) sum(v == max(v)))
    ok <- is.finite(best) & (ties == 1)
    assigned[ok] <- colnames(S)[winner[ok]]
    details$scores <- S

  } else if (strategy %in% c("partition_kmeans","partition_hclust")) {
    if (is.null(k_partition) || !is.finite(k_partition) || k_partition < 1L) {
      stop("`k_partition` must be a positive integer for partitioning strategies.")
    }
    ded <- .dedup_cols_by_key(Gcut)
    F <- ded$M
    if (!ncol(F)) stop("No informative features remain after deduplication.")
    F_sc <- if (isTRUE(feature_scale)) scale(F) else as.matrix(F)

    if (strategy == "partition_kmeans") {
      if (!is.null(kmeans_seed)) set.seed(as.integer(kmeans_seed))
      km <- stats::kmeans(F_sc, centers = as.integer(k_partition), iter.max = 100, nstart = 10)
      grp <- km$cluster
      assigned <- paste0("k", grp)
      names(assigned) <- rownames(F_sc)
      if (min_group_size > 1) {
        counts <- table(assigned)
        small  <- names(counts)[counts < min_group_size]
        if (length(small)) assigned[assigned %in% small] <- "not assigned"
      }
      details$partition <- data.frame(plot = names(assigned), cluster = assigned, stringsAsFactors = FALSE)

    } else {
      if (use_parallelDist && !requireNamespace("parallelDist", quietly = TRUE)) {
        warning("parallelDist not available; falling back to stats::dist().")
        use_parallelDist <- FALSE
      }
      d <- if (use_parallelDist) parallelDist::parDist(F_sc, method = dist_method)
      else stats::dist(F_sc, method = dist_method)

      if (use_fastcluster && !requireNamespace("fastcluster", quietly = TRUE)) {
        warning("fastcluster not available; falling back to stats::hclust().")
        use_fastcluster <- FALSE
      }
      hc <- if (use_fastcluster) fastcluster::hclust(d, method = clust_method)
      else stats::hclust(d, method = clust_method)

      grp_cut <- stats::cutree(hc, k = as.integer(k_partition))
      if (is.null(names(grp_cut))) names(grp_cut) <- rownames(F_sc)
      grp_cut <- grp_cut[rownames(F_sc)]
      grp_lbl <- .relabel_by_first_leaf(hc, grp_cut)
      assigned <- paste0("h", grp_lbl)
      names(assigned) <- rownames(F_sc)

      if (min_group_size > 1) {
        counts <- table(assigned)
        small  <- names(counts)[counts < min_group_size]
        if (length(small)) assigned[assigned %in% small] <- "not assigned"
      }

      details$partition <- data.frame(
        plot = names(grp_lbl),
        cluster = assigned,
        cluster_orig = grp_cut,
        stringsAsFactors = FALSE
      )
      details$hclust <- hc
      details$dist_method <- dist_method
      details$clust_method <- clust_method
    }
  } else {
    stop("Unknown strategy: ", strategy)
  }

  if (strategy %in% c("topology_only","topology_cover","cover_strict","cover_fuzzy")) {
    if (min_group_size > 1) {
      counts <- table(assigned)
      small  <- names(counts)[counts < min_group_size & names(counts) != "not assigned"]
      if (length(small)) assigned[assigned %in% small] <- "not assigned"
      details$collapsed_groups <- small
    } else {
      details$collapsed_groups <- character(0)
    }
  }

  structure(stats::setNames(assigned, names(assigned)), details = details)
}
