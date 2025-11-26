#' Group Cocktail clusters by fuzzy species membership
#'
#' @description
#' Given a set of Cocktail nodes (e.g. from \code{\link{clusters_at_cut}}),
#' recompute **direct (Option A)** \eqn{\phi} for all species versus the
#' **union groups** defined by these nodes (as in
#' \code{assign_releves}(strategy = "fuzzy")), then compute a dissimilarity
#' between groups based on their fuzzy species profiles and perform
#' hierarchical clustering.
#'
#' @param x A \code{"cocktail"} object (result of \code{\link{cocktail_cluster}}).
#' @param vegmatrix A matrix or data.frame with **plots in rows** and
#'   **species in columns**, containing the original vegetation data used for
#'   clustering. Row names must include (and will be reordered to) the plots
#'   in \code{x}; column names must include the species used in \code{x}.
#' @param clusters Cluster identifiers to be grouped. Can be:
#'   \itemize{
#'     \item a **vector** of node labels/IDs (e.g., \code{c("c_12","c_27")} or
#'       \code{c(12,27)}), each defining a separate group; or
#'     \item a **list**, where each element is a vector of node labels/IDs to be
#'       \strong{unioned} into a single group (OR of their plot memberships),
#'       e.g. \code{list(c("c_12","c_34"), 57)}.
#'   }
#' @param phi_fuzzy Numeric in (0,1); species with \eqn{\phi <} \code{phi_fuzzy}
#'   for a given group are set to 0 for that group. Defaults to 0.30. Use
#'   \code{phi_fuzzy = 0} to keep the full \eqn{\phi} vectors (no threshold).
#' @param min_species Integer ≥ 1. Minimum number of species in the
#'   **topological Cocktail group** (union of the selected nodes). Any group
#'   whose union has fewer than this many species is dropped \emph{before}
#'   computing fuzzy \eqn{\phi} and dissimilarities. Default 3.
#' @param dist_method Character; method for computing dissimilarities between
#'   groups based on their fuzzy species profiles. One of
#'   \code{c("bray","jaccard","euclidean","manhattan","correlation")}.
#'   \itemize{
#'     \item \code{"bray"}: Bray--Curtis dissimilarity.
#'     \item \code{"jaccard"}: 1 − Jaccard similarity (using fuzzy abundances).
#'     \item \code{"euclidean"}, \code{"manhattan"}: passed to \code{\link[stats]{dist}}.
#'     \item \code{"correlation"}: 1 − Pearson correlation between group profiles.
#'   }
#' @param hclust_method Character; linkage method passed to
#'   \code{\link[stats]{hclust}} (e.g. \code{"average"}, \code{"complete"},
#'   \code{"single"}, ...).
#'
#' @return
#' An \code{\link[stats]{hclust}} object. The result carries attributes:
#' \itemize{
#'   \item \code{"phi_matrix"} — the species × group \eqn{\phi} matrix
#'         (after applying \code{phi_fuzzy});
#'   \item \code{"groups"} — the final group labels used in clustering;
#'   \item \code{"clusters"} — the original \code{clusters} argument;
#'   \item \code{"species_per_group_topo"} — number of **topological**
#'         species per retained group (union of Cocktail nodes).
#' }
#'
#' @details
#' The function first builds a plots × groups membership matrix \eqn{G} by
#' taking the OR-union of the selected node columns in \code{x$Plot.cluster}
#' for each group, and a corresponding union of species membership from
#' \code{x$Cluster.species}. Groups whose topological species sets contain
#' fewer than \code{min_species} species are dropped.
#'
#' It then recomputes **direct (Option A)** \eqn{\phi} for all species
#' (from \code{vegmatrix}, binarized) against these group membership vectors,
#' exactly as in the `"fuzzy"` branch of \code{\link{assign_releves}}.
#' Species with \eqn{\phi <} \code{phi_fuzzy} for a group are set to 0; the
#' resulting columns (groups) are fuzzy species profiles. Dissimilarities
#' between groups are computed from these profiles using the chosen
#' \code{dist_method}, and hierarchical clustering is performed with
#' \code{\link[stats]{hclust}} and the requested \code{hclust_method}.
#'
#' Intuitively, \code{min_species} controls how “small” a Cocktail group
#' (in terms of its topological species set) is allowed to be, whereas
#' \code{phi_fuzzy} controls how strict you are about species–group fidelity
#' when constructing fuzzy profiles.
#'
#' @seealso \code{\link{clusters_at_cut}}, \code{\link{assign_releves}},
#'   \code{\link{cocktail_cluster}}
#' @importFrom stats dist hclust cor as.dist
#' @export
group_clusters_new <- function(
    x,
    vegmatrix,
    clusters,
    phi_fuzzy     = 0.30,
    min_species   = 3L,
    dist_method   = c("bray","jaccard","euclidean","manhattan","correlation"),
    hclust_method = "average"
) {
  ## ---- helpers ------------------------------------------------------------
  .is_cocktail <- function(obj) {
    is.list(obj) &&
      all(c("Cluster.species","Cluster.merged","Cluster.height","Plot.cluster") %in% names(obj))
  }

  # parse clusters arg (vector or list) into list of integer vectors
  .parse_clusters_arg <- function(v) {
    parse_one <- function(x) {
      if (is.character(x)) {
        as.integer(sub("^c_", "", x))
      } else {
        as.integer(x)
      }
    }
    if (is.list(v)) {
      lapply(v, parse_one)
    } else {
      lst <- as.list(v)
      lapply(lst, parse_one)
    }
  }

  # reduce node IDs to topmost within that set (drop descendants of any kept node)
  .keep_topmost_within <- function(ids, CM) {
    ids <- sort(unique(ids))
    if (!length(ids)) return(ids)
    kids <- unique(as.integer(CM[ids, , drop = FALSE][CM[ids, , drop = FALSE] > 0]))
    sort(setdiff(ids, intersect(ids, kids)))
  }

  # build groups from explicit clusters arg:
  # returns list with $group_nodes, $group_labels, $G (N×G membership),
  # and $n_species_topo (length G, number of topological species per group)
  .build_groups_from_clusters <- function(CS, CM, PC, clusters, min_species) {
    if (missing(clusters) || is.null(clusters)) {
      stop("`clusters` must be provided (e.g., output of clusters_at_cut()).")
    }
    node_groups <- .parse_clusters_arg(clusters)
    if (!length(node_groups)) stop("No valid cluster IDs found in `clusters`.")

    # validate & reduce each group to topmost nodes
    node_groups <- lapply(node_groups, function(ids) {
      ids <- ids[is.finite(ids) & ids > 0L]
      if (!length(ids)) return(integer(0))
      ids <- ids[ids <= nrow(CS)]
      .keep_topmost_within(ids, CM)
    })

    # discard groups that became empty
    non_empty <- vapply(node_groups, function(v) length(v) > 0L, logical(1))
    if (!any(non_empty)) stop("After ancestor/descendant reduction, no groups remain.")
    node_groups <- node_groups[non_empty]

    min_species <- as.integer(min_species)
    if (!is.finite(min_species) || min_species < 1L) min_species <- 1L

    Glist           <- list()
    lablist         <- character(0L)
    n_species_topo  <- integer(0L)

    for (g in seq_along(node_groups)) {
      ids <- node_groups[[g]]

      # topological species union for this group
      sp_union_logical <- colSums(CS[ids, , drop = FALSE] > 0L) > 0L
      n_sp <- sum(sp_union_logical)

      # skip groups with too few topological species
      if (n_sp < min_species) next

      cols <- paste0("c_", ids)
      missing_cols <- setdiff(cols, colnames(PC))
      if (length(missing_cols)) {
        stop("Internal mismatch: Plot.cluster is missing: ",
             paste(missing_cols, collapse = ", "))
      }
      Gi <- rowSums(PC[, cols, drop = FALSE] > 0) > 0

      Glist[[length(Glist) + 1L]] <- as.numeric(Gi)
      lablist[length(lablist) + 1L] <- paste0("g_", paste(ids, collapse = "_"))
      n_species_topo[length(n_species_topo) + 1L] <- n_sp
    }

    if (!length(Glist)) {
      stop("After applying min_species on topological species sets, no groups remain.")
    }

    G <- do.call(cbind, Glist)
    rownames(G) <- rownames(PC)
    colnames(G) <- lablist

    list(
      group_nodes      = node_groups,
      group_labels     = lablist,
      G                = G,
      n_species_topo   = n_species_topo
    )
  }

  # Direct (Option A) phi for species vs groups given X (N×S 0/1) and G (N×G 0/1)
  .phi_direct_species_by_groups <- function(X, G) {
    storage.mode(X) <- "double"
    storage.mode(G) <- "double"
    N  <- nrow(X); S <- ncol(X); Gk <- ncol(G)
    a  <- crossprod(X, G)                  # S×G
    p  <- colSums(X)                       # length S
    g1 <- colSums(G)                       # length G
    b  <- matrix(p,  nrow = S, ncol = Gk) - a
    c  <- matrix(g1, nrow = S, ncol = Gk, byrow = TRUE) - a
    d  <- (N - matrix(g1, nrow = S, ncol = Gk, byrow = TRUE)) - b
    den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
    Phi <- (a * d - b * c) / den
    Phi[!is.finite(den) | den <= 0] <- 0
    Phi
  }

  # pairwise dissimilarity between rows of M (groups × species)
  .group_dissimilarity <- function(M, method) {
    if (method %in% c("euclidean","manhattan")) {
      return(stats::dist(M, method = method))
    }
    if (method == "correlation") {
      cc <- stats::cor(t(M), use = "pairwise.complete.obs")
      return(stats::as.dist(1 - cc))
    }

    # custom Bray–Curtis and Jaccard for fuzzy data
    n <- nrow(M)
    if (n < 2L) stop("Need at least two groups to compute dissimilarities.")
    D <- matrix(0, n, n)
    for (i in 1L:(n - 1L)) {
      xi <- M[i, ]
      for (j in (i + 1L):n) {
        xj <- M[j, ]
        if (method == "bray") {
          num <- sum(abs(xi - xj))
          den <- sum(xi + xj)
          d   <- if (den > 0) num / den else 0
        } else { # jaccard
          num <- sum(pmin(xi, xj))
          den <- sum(pmax(xi, xj))
          d   <- if (den > 0) 1 - num / den else 0
        }
        D[i, j] <- D[j, i] <- d
      }
    }
    dimnames(D) <- list(rownames(M), rownames(M))
    stats::as.dist(D)
  }

  ## ---- checks / setup -----------------------------------------------------
  if (!.is_cocktail(x)) {
    stop("`x` must be a Cocktail object with components ",
         "Cluster.species, Cluster.merged, Cluster.height, and Plot.cluster.")
  }

  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  PC <- x$Plot.cluster

  # normalize Plot.cluster colnames like in assign_releves()
  if (is.null(colnames(PC))) {
    colnames(PC) <- paste0("c_", seq_len(ncol(PC)))
  } else {
    num_only <- grepl("^\\d+$", colnames(PC))
    colnames(PC)[num_only] <- paste0("c_", colnames(PC)[num_only])
  }

  # build groups from clusters arg (unions of nodes) and filter by min_species (topological)
  grp           <- .build_groups_from_clusters(CS, CM, PC, clusters, min_species)
  group_labels  <- grp$group_labels
  Ggrp          <- as.matrix(grp$G)   # N × G
  n_sp_topo     <- grp$n_species_topo

  # vegmatrix: align plots and species
  if (is.null(vegmatrix)) {
    stop("`vegmatrix` is required to compute species × group phi.")
  }
  vm <- as.matrix(vegmatrix)
  vm[is.na(vm)] <- 0

  plot_names <- rownames(PC)
  if (is.null(plot_names)) {
    stop("`x$Plot.cluster` must have row names (plot IDs).")
  }
  if (is.null(rownames(vm))) {
    stop("`vegmatrix` must have row names matching the plot IDs in `x`.")
  }
  if (!all(plot_names %in% rownames(vm))) {
    stop("Row names of `vegmatrix` must include all plots in `x`.")
  }
  vm <- vm[plot_names, , drop = FALSE]

  spp_all <- colnames(CS)
  if (is.null(spp_all)) {
    stop("`x$Cluster.species` must have species column names.")
  }
  common <- intersect(colnames(vm), spp_all)
  if (!length(common)) {
    stop("No overlapping species between `vegmatrix` and Cocktail species.")
  }
  vm <- vm[, common, drop = FALSE]

  # presence/absence matrix (N × S')
  X <- (vm > 0)
  storage.mode(X) <- "double"

  # recompute direct phi (species × group)
  Phi <- .phi_direct_species_by_groups(X, Ggrp)
  rownames(Phi) <- colnames(X)
  colnames(Phi) <- group_labels

  # threshold by phi_fuzzy (phi_fuzzy = 0 means no threshold)
  phi_fuzzy <- max(0, min(1, as.numeric(phi_fuzzy)))
  if (phi_fuzzy > 0) {
    Phi[Phi < phi_fuzzy] <- 0
  }

  # groups × species matrix for distances
  M <- t(Phi)
  if (nrow(M) < 2L) {
    stop("After applying min_species on topological species sets, fewer than two groups remain; ",
         "cannot build an hclust object.")
  }

  dist_method <- match.arg(dist_method)
  D <- .group_dissimilarity(M, dist_method)

  ## ---- hclust -------------------------------------------------------------
  hc <- stats::hclust(D, method = hclust_method)
  attr(hc, "phi_matrix")           <- Phi
  attr(hc, "groups")               <- group_labels
  attr(hc, "clusters")             <- clusters
  attr(hc, "species_per_group_topo") <- n_sp_topo
  hc
}
