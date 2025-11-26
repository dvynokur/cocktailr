#' Assign relevés (plots) to Cocktail groups using relative cover and species–cluster φ
#'
#' @description
#' Assign plots to groups derived from a Cocktail tree using only the
#' information stored in a fitted Cocktail object:
#' \code{x$Plot.cluster} (relative covers per plot × node),
#' \code{x$Cluster.species} (topological species per node),
#' and optionally \code{x$Species.cluster.phi} (species–cluster φ).
#'
#' Four strategies are supported:
#'
#' - **"num_sp"** — For each group, count its *topological* species
#'   (from \code{x$Cluster.species}). For any plot, the eligible groups
#'   are those to which the plot belongs (non-zero membership in
#'   \code{x$Plot.cluster}, or OR-union for union groups). The winner is
#'   the group with the **largest number of topological species** among
#'   those eligible (ties → `"not assigned"`).
#'
#' - **"cover"** — For each group, the score for a plot is the
#'   **sum of relative covers** from \code{x$Plot.cluster} across all
#'   nodes in that group (for single-node groups this is just the
#'   relative cover of that node). The winner is the group with
#'   **highest relative cover** for that plot (ties → `"not assigned"`).
#'
#' - **"phi_topo"** — Requires \code{x$Species.cluster.phi}. For each
#'   group, build a **topological species set** (union of its nodes'
#'   topological species) and a **species–group φ** by taking, for each
#'   species, the **maximum φ** across the group's nodes (negative φ are
#'   set to 0). The group-level “φ-strength” is the **sum of these φ**
#'   over its topological species. For a plot, eligible groups are those
#'   to which the plot belongs; the winner has the **largest φ-strength**
#'   among eligible groups.
#'
#' - **"phi_cover"** — Requires \code{x$Species.cluster.phi}. As in
#'   `"phi_topo"`, compute species–group φ (max over nodes in group),
#'   but then **threshold** by \code{min_phi}: φ < \code{min_phi} are
#'   set to 0. For each group, define its φ-strength as the sum of the
#'   remaining φ over its topological species. The score for a plot and
#'   group is then **relative cover × φ-strength**. The winner is the
#'   group with the highest such score among eligible groups.
#'
#' In all strategies, a plot is considered **eligible** for a group if
#' it has any non-zero membership for the group's nodes in
#' \code{x$Plot.cluster} (OR-union for union groups). Per-plot scores
#' are compared only across eligible groups. Ties or all-NA scores yield
#' `"not assigned"`.
#'
#' @param x A Cocktail object (list) with components
#'   \code{Cluster.species}, \code{Cluster.merged},
#'   \code{Cluster.height}, and \code{Plot.cluster}, e.g. from
#'   \code{\link{cocktail_cluster}()}. For strategies
#'   \code{"phi_topo"} and \code{"phi_cover"}, the object must also
#'   contain \code{Species.cluster.phi}, a matrix of Option A φ values
#'   (species × nodes).
#' @param strategy Character, one of
#'   \code{c("num_sp","cover","phi_topo","phi_cover")}.
#' @param phi_cut Numeric in (0,1). Height cut used to choose **parent
#'   clusters** (unless \code{clusters} is supplied). Parent clusters at
#'   a cut are the topmost nodes with height ≥ \code{phi_cut} that do
#'   not have an ancestor also at/above the cut.
#' @param clusters Optional selection of clusters **instead of**
#'   \code{phi_cut}.
#'   \itemize{
#'     \item If a **vector** of node labels/IDs
#'       (e.g., \code{c("c_123","c_456")} or \code{c(123,456)}),
#'       each element forms a separate group.
#'     \item If a **list**, each element is a **union group** of several
#'       nodes, e.g., \code{list(c("c_123","c_234"), 567)}.
#'   }
#'   Within each group, if both an ancestor and its descendant are
#'   present, only the **topmost** (ancestor) is kept to avoid
#'   double-counting.
#' @param min_phi Numeric in (0,1). Minimum φ value for
#'   \code{strategy = "phi_cover"}: species–group φ < \code{min_phi} are
#'   set to 0 before computing group φ-strength. Default 0.2. Ignored
#'   for other strategies.
#' @param min_group_size Integer ≥ 1. After assignment, any group with
#'   fewer than this many plots is relabeled `"not assigned"`. Default 1
#'   (no collapsing).
#'
#' @details
#' \strong{Groups and union groups.}
#' Groups are defined either by \code{phi_cut} (parent nodes at that cut)
#' or explicitly by \code{clusters}. When \code{clusters} is a list, each
#' element is a union group of several nodes. Plot membership in a union
#' group is taken as the OR-union of membership over its nodes
#' (non-zero entries in \code{x$Plot.cluster}); the topological species
#' of a union group are the union of the nodes' topological species
#' (rows of \code{x$Cluster.species}).
#'
#' \strong{Species–group φ for union groups.}
#' For \code{"phi_topo"} and \code{"phi_cover"}, species–group φ for a
#' union group is derived from \code{x$Species.cluster.phi} by taking the
#' maximum φ over the group's nodes for each species. Negative φ values
#' are set to 0 before aggregation. For \code{"phi_cover"}, φ-values
#' below \code{min_phi} are additionally set to 0.
#'
#' @return
#' A named character vector of length equal to the number of plots, with
#' values equal to group labels (e.g. \code{"g_12"}, \code{"g_5_7"}) or
#' \code{"not assigned"}. An attribute \code{"details"} is attached
#' containing the strategy, groups used, scoring summaries, and the
#' arguments used.
#'
#' @seealso \code{\link{cocktail_cluster}}, \code{\link{clusters_at_cut}}
#' @export
assign_releves <- function(
    x,
    strategy       = c("num_sp","cover","phi_topo","phi_cover"),
    phi_cut        = 0.30,
    clusters       = NULL,   # vector or list (for unions)
    min_phi        = 0.2,
    min_group_size = 1L
) {
  ## ---- helpers ------------------------------------------------------------
  .is_cocktail <- function(obj) {
    is.list(obj) &&
      all(c("Cluster.species","Cluster.merged","Cluster.height","Plot.cluster") %in% names(obj))
  }

  # parents (topmost) at a cut
  .parents_at_cut <- function(CM, H, cut) {
    idx <- which(H >= cut)
    if (!length(idx)) return(integer(0))
    kids <- CM[idx, , drop = FALSE]
    children <- unique(as.integer(kids[kids > 0]))
    sort(setdiff(idx, intersect(idx, children)))
  }

  # parse clusters arg (vector or list); returns a list of integer vectors (each = one group)
  .parse_clusters_arg <- function(v) {
    if (is.null(v)) return(NULL)
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

  # reduce a set of node IDs to topmost within that set (drop descendants of any kept node)
  .keep_topmost_within <- function(ids, CM) {
    ids <- sort(unique(ids))
    if (!length(ids)) return(ids)
    kids <- unique(as.integer(CM[ids, , drop = FALSE][CM[ids, , drop = FALSE] > 0]))
    sort(setdiff(ids, intersect(ids, kids)))
  }

  # Build groups: from phi_cut or clusters arg.
  # Returns:
  #   group_nodes      — list of integer vectors (node IDs per group)
  #   group_labels     — character vector of group labels ("g_<id>" or "g_id1_id2")
  #   G                — N × G membership (0/1)
  #   cover_mat        — N × G relative cover (sum over nodes in group)
  #   species_sets_topo — list of character vectors (topological species per group)
  .build_groups <- function(CS, CM, H, PC, phi_cut, clusters) {
    spp_all <- colnames(CS)
    N <- nrow(PC)

    # choose node sets per group
    node_groups <- if (is.null(clusters)) {
      as.list(.parents_at_cut(CM, H, phi_cut))
    } else {
      .parse_clusters_arg(clusters)
    }

    if (!length(node_groups)) {
      stop("No groups available (check phi_cut or clusters).")
    }

    # validate & reduce each group to topmost
    node_groups <- lapply(seq_along(node_groups), function(i) {
      ids <- node_groups[[i]]
      if (any(!is.finite(ids))) stop("`clusters` contains non-integer/invalid IDs.")
      ids <- ids[ids > 0]
      .keep_topmost_within(ids, CM)
    })

    # discard empty groups (after reduction)
    keep_group <- vapply(node_groups, function(v) length(v) > 0, logical(1))
    if (!any(keep_group)) {
      stop("After ancestor/descendant reduction, no groups remain.")
    }
    node_groups <- node_groups[keep_group]

    # membership and cover per group: OR and sum over Plot.cluster columns
    Glist           <- list()
    Covlist         <- list()
    lablist         <- character(length(node_groups))
    species_sets    <- vector("list", length(node_groups))

    for (g in seq_along(node_groups)) {
      ids  <- node_groups[[g]]
      cols <- paste0("c_", ids)
      missing_cols <- setdiff(cols, colnames(PC))
      if (length(missing_cols)) {
        stop("Internal mismatch: Plot.cluster is missing: ",
             paste(missing_cols, collapse = ", "))
      }

      # subset (possibly sparse) and convert once
      sub_pc <- as.matrix(PC[, cols, drop = FALSE])
      # membership: any non-zero entry across nodes
      Gi  <- rowSums(sub_pc != 0) > 0
      Glist[[g]] <- as.numeric(Gi)
      # cover: sum of relative cover across nodes in the group
      Covlist[[g]] <- rowSums(sub_pc)

      lablist[g] <- if (length(ids) == 1L) {
        paste0("g_", ids)
      } else {
        paste0("g_", paste(ids, collapse = "_"))
      }

      # topological species union across these node IDs
      sp_union_logical <- apply(CS[ids, , drop = FALSE] > 0, 2, any)
      species_sets[[g]] <- spp_all[sp_union_logical]
    }

    G   <- do.call(cbind, Glist)
    Cov <- do.call(cbind, Covlist)

    colnames(G)   <- lablist
    colnames(Cov) <- lablist

    list(
      group_nodes       = node_groups,
      group_labels      = lablist,
      G                 = G,
      cover_mat         = Cov,
      species_sets_topo = species_sets
    )
  }

  ## ---- checks / setup -----------------------------------------------------
  strategy <- match.arg(strategy)

  if (!.is_cocktail(x)) {
    stop("`x` must be a Cocktail object with Cluster.species, Cluster.merged, Cluster.height, and Plot.cluster.")
  }

  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  PC <- x$Plot.cluster

  # normalise Plot.cluster colnames to "c_<index>"
  if (is.null(colnames(PC))) {
    colnames(PC) <- paste0("c_", seq_len(ncol(PC)))
  } else {
    num_only <- grepl("^\\d+$", colnames(PC))
    colnames(PC)[num_only] <- paste0("c_", colnames(PC)[num_only])
  }

  plot_names <- rownames(PC)
  if (is.null(plot_names)) {
    plot_names <- paste0("plot_", seq_len(nrow(PC)))
    rownames(PC) <- plot_names
  }

  spp_all <- colnames(CS)
  if (is.null(spp_all)) {
    stop("`x$Cluster.species` must have species column names.")
  }

  # build groups (nodes, labels, membership, cover, topological species sets)
  grp <- .build_groups(CS, CM, H, PC, phi_cut, clusters)
  group_labels      <- grp$group_labels
  Ggrp              <- as.matrix(grp$G)           # N × G, 0/1
  cover_mat         <- as.matrix(grp$cover_mat)   # N × G, numeric
  species_sets_topo <- grp$species_sets_topo     # list of length G

  n_plots <- nrow(Ggrp)
  n_groups <- ncol(Ggrp)

  # hygiene
  min_group_size <- as.integer(min_group_size)
  if (!is.finite(min_group_size) || min_group_size < 1L) {
    min_group_size <- 1L
  }
  min_phi <- max(0, min(1, as.numeric(min_phi)))

  # base output
  assigned <- rep("not assigned", n_plots)
  names(assigned) <- plot_names

  details <- list(
    strategy       = strategy,
    phi_cut        = phi_cut,
    clusters       = clusters,
    groups_used    = group_labels,
    min_phi        = min_phi,
    min_group_size = min_group_size
  )

  ## ---- precompute φ-based summaries if needed -----------------------------
  Phi_pos <- NULL
  group_phi_topo  <- NULL
  group_phi_cover <- NULL

  if (strategy %in% c("phi_topo","phi_cover")) {
    if (!("Species.cluster.phi" %in% names(x)) || is.null(x$Species.cluster.phi)) {
      stop(
        "`x$Species.cluster.phi` is not available.\n",
        "Recompute the Cocktail clustering with `species_cluster_phi = TRUE`, e.g.:\n",
        "  x <- cocktail_cluster(vegmatrix, species_cluster_phi = TRUE, ...)\n",
        "and then call `assign_releves()` with strategy = 'phi_topo' or 'phi_cover'."
      )
    }
    Phi <- x$Species.cluster.phi
    if (!is.matrix(Phi)) stop("`x$Species.cluster.phi` must be a matrix.")

    if (is.null(rownames(Phi))) {
      stop("`x$Species.cluster.phi` must have species row names.")
    }

    # align species
    if (!all(spp_all %in% rownames(Phi))) {
      missing_sp <- setdiff(spp_all, rownames(Phi))
      stop("`Species.cluster.phi` is missing species: ",
           paste(head(missing_sp, 10), collapse = ", "),
           if (length(missing_sp) > 10) " ..." else "")
    }
    Phi <- Phi[spp_all, , drop = FALSE]  # reorder rows
    Phi_pos <- pmax(Phi, 0)
    storage.mode(Phi_pos) <- "double"

    node_colnames <- colnames(Phi_pos)
    if (is.null(node_colnames)) {
      stop("`Species.cluster.phi` must have column names (node IDs, e.g. 'c_1').")
    }

    # helper: map node ID to column index in Phi_pos
    .node_to_col_idx <- function(k) {
      cand1 <- paste0("c_", k)
      if (cand1 %in% node_colnames) return(match(cand1, node_colnames))
      cand2 <- as.character(k)
      if (cand2 %in% node_colnames) return(match(cand2, node_colnames))
      NA_integer_
    }

    group_phi_topo  <- numeric(n_groups)
    group_phi_cover <- numeric(n_groups)
    names(group_phi_topo)  <- group_labels
    names(group_phi_cover) <- group_labels

    for (g in seq_len(n_groups)) {
      ids <- grp$group_nodes[[g]]
      if (!length(ids)) next

      col_idx <- vapply(ids, .node_to_col_idx, integer(1L))
      if (anyNA(col_idx)) {
        missing_ids <- ids[is.na(col_idx)]
        stop("`Species.cluster.phi` is missing columns for nodes: ",
             paste(missing_ids, collapse = ", "))
      }

      sp_set <- species_sets_topo[[g]]
      if (!length(sp_set)) {
        group_phi_topo[g]  <- 0
        group_phi_cover[g] <- 0
        next
      }

      Phi_g <- Phi_pos[sp_set, col_idx, drop = FALSE]  # |S_g| × |ids|
      # species–group φ = max over nodes in group
      phi_s <- apply(Phi_g, 1L, max)

      # for phi_topo: just sum of non-negative φ
      group_phi_topo[g] <- sum(phi_s, na.rm = TRUE)

      # for phi_cover: apply min_phi threshold first
      phi_s_thr <- phi_s
      phi_s_thr[phi_s_thr < min_phi] <- 0
      group_phi_cover[g] <- sum(phi_s_thr, na.rm = TRUE)
    }

    details$group_phi_topo  <- group_phi_topo
    details$group_phi_cover <- group_phi_cover
  }

  ## ---- compute scores per strategy ----------------------------------------
  # score matrix S: N × G
  S <- matrix(-Inf, nrow = n_plots, ncol = n_groups,
              dimnames = list(plot_names, group_labels))

  if (strategy == "num_sp") {
    # number of topological species per group
    sp_counts <- vapply(species_sets_topo, length, integer(1L))
    for (g in seq_len(n_groups)) {
      eligible <- Ggrp[, g] > 0
      S[eligible, g] <- sp_counts[g]
    }

  } else if (strategy == "cover") {
    for (g in seq_len(n_groups)) {
      eligible <- Ggrp[, g] > 0
      S[eligible, g] <- cover_mat[eligible, g]
    }

  } else if (strategy == "phi_topo") {
    for (g in seq_len(n_groups)) {
      eligible <- Ggrp[, g] > 0
      S[eligible, g] <- group_phi_topo[g]
    }

  } else if (strategy == "phi_cover") {
    for (g in seq_len(n_groups)) {
      eligible <- Ggrp[, g] > 0 & cover_mat[, g] > 0
      S[eligible, g] <- cover_mat[eligible, g] * group_phi_cover[g]
    }
  }

  ## ---- pick winners per plot ----------------------------------------------
  # max.col will always return a column, but we guard using 'best' and 'ties'
  winner <- max.col(S, ties.method = "first")
  best   <- S[cbind(seq_len(nrow(S)), winner)]
  ties   <- apply(S, 1L, function(v) sum(v == max(v)))

  ok <- is.finite(best) & (ties == 1L)
  assigned[ok] <- group_labels[winner[ok]]

  ## ---- collapse small groups ----------------------------------------------
  if (min_group_size > 1L) {
    counts <- table(assigned)
    small  <- names(counts)[counts < min_group_size & names(counts) != "not assigned"]
    if (length(small)) {
      assigned[assigned %in% small] <- "not assigned"
    }
    details$collapsed_groups <- small
  } else {
    details$collapsed_groups <- character(0)
  }

  structure(stats::setNames(assigned, names(assigned)), details = details)
}
