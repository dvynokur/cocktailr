#' Assign releves (plots) to Cocktail groups using covers and species-cluster phi
#'
#' @description
#' Assign plots to groups derived from a Cocktail tree using:
#' \itemize{
#'   \item \code{x$Plot.cluster} - cluster membership per plot (non-zero = member);
#'   \item \code{x$Cluster.species} - node-constituent species per node;
#'   \item \code{vegmatrix} - per-plot species covers (plots x species);
#'   \item optionally \code{x$Species.cluster.phi} - species-cluster \eqn{\phi}.
#' }
#'
#' Strategies (per plot, per group):
#'
#' \describe{
#'
#' \item{\code{"count"}}{
#'   For each cluster, count how many of its node-constituent species are present
#'   in the releve. The winner is the cluster with the largest count.
#'   If several clusters share the same maximum, the plot is labeled \code{"+"}.
#'   If no cluster has any species present, the plot gets \code{NA}.
#' }
#'
#' \item{\code{"cover"}}{
#'   For each cluster, sum the covers of its node-constituent species that are
#'   present in the releve. The winner is the cluster with the highest sum.
#'   If there is a tie, use the number of species (as in \code{"count"})
#'   as a tie-breaker; if still tied, label \code{"+"}. Plots with zero
#'   cover for all clusters get \code{NA}.
#'
#'   If \code{x$Plot.cluster} appears to be strictly binary (0/1), the function
#'   issues a warning and automatically falls back to \code{strategy = "count"}.
#' }
#'
#' \item{\code{"phi_node"}}{
#'   Requires \code{x$Species.cluster.phi}. For each cluster and species,
#'   define \eqn{\phi(s,g)} as the maximum of \eqn{\phi} across that group's
#'   nodes (negatives set to 0). For a given plot, restrict to the cluster's
#'   node-constituent species that are present in the releve and sum their
#'   \eqn{\phi} values. The winner has the largest such sum. Ties are broken
#'   first by total cover of these species, then by their species count;
#'   if still tied, label \code{"+"}. If all sums are 0, the plot gets \code{NA}.
#' }
#'
#' \item{\code{"phi_cover_node"}}{
#'   Requires \code{x$Species.cluster.phi}. As in \code{"phi_node"}, use
#'   the cluster's node-constituent species only, but score each cluster in a plot
#'   by the sum of \eqn{\mathrm{cover}(i,s)\,\phi(s,g)} over its node-constituent
#'   species present in the releve. The winner has the largest score.
#'   If there is a tie, use the number of these species as a tie-breaker;
#'   if still tied, label \code{"+"}. Plots with all scores 0 get \code{NA}.
#' }
#'
#' \item{\code{"phi_cover"}}{
#'   Requires \code{x$Species.cluster.phi}. For each cluster and species,
#'   define \eqn{\phi(s,g)} as above, then keep only species with
#'   \eqn{\phi(s,g) \ge \code{min_phi}}. For a plot, restrict to these
#'   species that are present in the releve and compute the score:
#'   \eqn{\sum \mathrm{cover}(i,s)\,\phi(s,g)}. The winner has the largest
#'   score; ties are broken by total cover of these species, then by their
#'   species count; if still tied, label \code{"+"}. If all scores are 0,
#'   the plot gets \code{NA}.
#' }
#'
#' \item{\code{"phi"}}{
#'   Requires \code{x$Species.cluster.phi}. For each cluster and species,
#'   define \eqn{\phi(s,g)} and keep only species with
#'   \eqn{\phi(s,g) \ge \code{min_phi}}. For a plot, restrict to these species
#'   that are present in the releve and sum their \eqn{\phi} values.
#'   The winner has the largest sum. Ties are broken by total cover of these
#'   species, then by their species count; if still tied, label \code{"+"}.
#'   If all sums are 0, the plot gets \code{NA}.
#' }
#'
#' }
#'
#' If a phi-based strategy is requested but \code{x$Species.cluster.phi} is
#' missing, the function issues a warning and automatically falls back to
#' \code{"count"} (when \code{x$Plot.cluster} is binary) or \code{"cover"}
#' (when it is not).
#'
#' Groups are defined either by \code{phi_cut} (parent nodes at that cut)
#' or explicitly by \code{clusters}, which can define union groups of nodes.
#'
#' @param x A Cocktail object (list) with components
#'   \code{Cluster.species}, \code{Cluster.merged}, \code{Cluster.height},
#'   and \code{Plot.cluster} (e.g. from \code{\link{cocktail_cluster}()}).
#'   For strategies \code{"phi_node"}, \code{"phi_cover_node"},
#'   \code{"phi_cover"}, and \code{"phi"}, \code{x$Species.cluster.phi}
#'   should be present; otherwise a fallback is used.
#' @param vegmatrix Numeric matrix or data.frame of covers (plots x species).
#'   Row names must include the plots in \code{x}; column names must include
#'   the species used in \code{x}. Missing values are treated as 0.
#' @param strategy One of
#'   \code{c("count","cover","phi_node","phi_cover_node","phi_cover","phi")}.
#' @param phi_cut Numeric in (0,1). Height cut used to choose parent clusters
#'   (topmost nodes with height >= \code{phi_cut}) when \code{clusters} is not given.
#' @param clusters Optional selection of clusters instead of \code{phi_cut}.
#'   If a vector of labels/IDs (e.g. \code{c("c_12","c_27")} or \code{c(12,27)}),
#'   each element defines a separate group. If a list, each element is a union
#'   of several nodes (OR of their memberships). Within each group, if both
#'   an ancestor and descendant are present, only the topmost (ancestor) is kept.
#' @param min_phi Numeric in (0,1). Minimum \eqn{\phi} for strategies
#'   \code{"phi_cover"} and \code{"phi"}: species-group \eqn{\phi < \code{min_phi}}
#'   are set to 0. Ignored by other strategies. Default 0.2.
#' @param min_group_size Integer >= 1. After assignment, any group with fewer
#'   than this many plots is relabeled \code{"-"}. Default 1.
#'
#' @return
#' A named character vector, with names = plot IDs and values = group labels
#' (e.g. \code{"g_12"}, \code{"g_5_7"}), \code{"+"} for ties,
#' \code{"-"} for groups collapsed by \code{min_group_size}, or \code{NA}
#' when no group wins. An attribute \code{"details"} is attached with
#' diagnostics (groups used, strategy, counts, etc.).
#'
#' @export
assign_releves <- function(
    x,
    vegmatrix,
    strategy       = c("count","cover","phi_node","phi_cover_node","phi_cover","phi"),
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

  .parents_at_cut <- function(CM, H, cut) {
    idx <- which(H >= cut)
    if (!length(idx)) return(integer(0))
    kids <- CM[idx, , drop = FALSE]
    children <- unique(as.integer(kids[kids > 0]))
    sort(setdiff(idx, intersect(idx, children)))
  }

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

  .keep_topmost_within <- function(ids, CM) {
    ids <- sort(unique(ids))
    if (!length(ids)) return(ids)
    kids <- unique(as.integer(CM[ids, , drop = FALSE][CM[ids, , drop = FALSE] > 0]))
    sort(setdiff(ids, intersect(ids, kids)))
  }

  .build_groups <- function(CS, CM, H, PC, phi_cut, clusters) {
    spp_all <- colnames(CS)

    node_groups <- if (is.null(clusters)) {
      as.list(.parents_at_cut(CM, H, phi_cut))
    } else {
      .parse_clusters_arg(clusters)
    }

    if (!length(node_groups)) {
      stop("No groups available (check phi_cut or clusters).")
    }

    node_groups <- lapply(node_groups, function(ids) {
      if (any(!is.finite(ids))) stop("`clusters` contains non-integer/invalid IDs.")
      ids <- ids[ids > 0]
      .keep_topmost_within(ids, CM)
    })

    keep_group <- vapply(node_groups, function(v) length(v) > 0, logical(1))
    if (!any(keep_group)) {
      stop("After ancestor/descendant reduction, no groups remain.")
    }
    node_groups <- node_groups[keep_group]

    Glist        <- list()
    lablist      <- character(length(node_groups))
    species_sets <- vector("list", length(node_groups))

    for (g in seq_along(node_groups)) {
      ids  <- node_groups[[g]]
      cols <- paste0("c_", ids)
      missing_cols <- setdiff(cols, colnames(PC))
      if (length(missing_cols)) {
        stop("Internal mismatch: Plot.cluster is missing: ",
             paste(missing_cols, collapse = ", "))
      }
      sub_pc <- as.matrix(PC[, cols, drop = FALSE])
      Gi <- rowSums(sub_pc != 0) > 0
      Glist[[g]] <- as.numeric(Gi)

      lablist[g] <- if (length(ids) == 1L) {
        paste0("g_", ids)
      } else {
        paste0("g_", paste(ids, collapse = "_"))
      }

      sp_union_logical <- apply(CS[ids, , drop = FALSE] > 0, 2, any)
      species_sets[[g]] <- spp_all[sp_union_logical]
    }

    G <- do.call(cbind, Glist)
    colnames(G) <- lablist

    list(
      group_nodes       = node_groups,
      group_labels      = lablist,
      G                 = G,
      species_sets_topo = species_sets
    )
  }

  .is_binary_cluster <- function(PC) {
    if (inherits(PC, "Matrix")) {
      nz <- PC@x
      if (!length(nz)) {
        TRUE
      } else {
        all(nz %in% c(0, 1))
      }
    } else {
      vals <- unique(as.vector(PC))
      vals <- vals[!is.na(vals)]
      if (!length(vals)) TRUE else all(vals %in% c(0, 1))
    }
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

  # check if Plot.cluster is effectively binary (0/1)
  is_binary_PC <- .is_binary_cluster(PC)

  # phi-based strategies: if phi matrix missing -> warn and fallback
  if (strategy %in% c("phi_node","phi_cover_node","phi_cover","phi") &&
      (!("Species.cluster.phi" %in% names(x)) || is.null(x$Species.cluster.phi))) {

    fallback <- if (is_binary_PC) "count" else "cover"
    warning(
      "Strategy '", strategy, "' requires x$Species.cluster.phi, which is not present.\n",
      "Recompute cocktail_cluster(..., species_cluster_phi = TRUE) to enable phi-based assignment.\n",
      "Falling back to strategy '", fallback, "'."
    )
    strategy <- fallback
  }

  # cover strategy with binary Plot.cluster -> warn and fallback to count
  if (strategy == "cover" && is_binary_PC) {
    warning(
      "Strategy 'cover' requested, but x$Plot.cluster appears to be binary (0/1).\n",
      "If cover values are available, recompute cocktail_cluster(..., plot_values = \"rel_cover\") ",
      "to enable cover-based assignment. Falling back to strategy 'count'."
    )
    strategy <- "count"
  }

  spp_all <- colnames(CS)
  if (is.null(spp_all)) {
    stop("`x$Cluster.species` must have species column names.")
  }

  if (missing(vegmatrix) || is.null(vegmatrix)) {
    stop("`vegmatrix` (covers) must be provided.")
  }
  vm <- as.matrix(vegmatrix)
  vm[is.na(vm)] <- 0
  storage.mode(vm) <- "double"

  if (is.null(colnames(vm))) {
    stop("`vegmatrix` must have species (column names).")
  }
  if (is.null(rownames(vm))) {
    stop("`vegmatrix` must have plot (row names) matching `x$Plot.cluster`.")
  }

  if (!all(plot_names %in% rownames(vm))) {
    missing_pl <- setdiff(plot_names, rownames(vm))
    stop("`vegmatrix` is missing plots used in Cocktail: ",
         paste(head(missing_pl, 10), collapse = ", "),
         if (length(missing_pl) > 10) " ..." else "")
  }
  vm <- vm[plot_names, , drop = FALSE]

  common <- intersect(colnames(vm), spp_all)
  if (!length(common)) {
    stop("No overlapping species between `vegmatrix` and Cocktail species.")
  }
  vm <- vm[, common, drop = FALSE]
  spp_all <- common
  CS <- CS[, spp_all, drop = FALSE]

  N <- nrow(vm)

  grp <- .build_groups(CS, CM, H, PC, phi_cut, clusters)
  group_labels      <- grp$group_labels
  Ggrp              <- as.matrix(grp$G)           # N x G, 0/1
  species_sets_topo <- grp$species_sets_topo     # list of length G

  n_groups <- ncol(Ggrp)

  min_group_size <- as.integer(min_group_size)
  if (!is.finite(min_group_size) || min_group_size < 1L) {
    min_group_size <- 1L
  }
  min_phi <- max(0, min(1, as.numeric(min_phi)))

  assigned <- rep(NA_character_, N)
  names(assigned) <- plot_names

  species_names <- colnames(vm)
  present_list <- lapply(seq_len(N), function(i) {
    species_names[ vm[i, ] > 0 ]
  })

  ## ---- phi precomputation if needed -----------------------------------------
  Phi_group_list <- NULL
  Phi_sets_phi   <- NULL

  if (strategy %in% c("phi_node","phi_cover_node","phi_cover","phi")) {
    Phi <- x$Species.cluster.phi
    if (is.null(Phi) || !is.matrix(Phi)) {
      stop("Internal error: phi-based strategy chosen but `x$Species.cluster.phi` is missing or not a matrix.")
    }

    if (is.null(rownames(Phi))) {
      stop("`x$Species.cluster.phi` must have species row names.")
    }

    if (!all(spp_all %in% rownames(Phi))) {
      missing_sp <- setdiff(spp_all, rownames(Phi))
      stop("`Species.cluster.phi` is missing species: ",
           paste(head(missing_sp, 10), collapse = ", "),
           if (length(missing_sp) > 10) " ..." else "")
    }
    Phi <- Phi[spp_all, , drop = FALSE]

    node_colnames <- colnames(Phi)
    if (is.null(node_colnames)) {
      stop("`Species.cluster.phi` must have column names (node IDs, e.g. 'c_1').")
    }

    .node_to_col_idx <- function(k) {
      cand1 <- paste0("c_", k)
      if (cand1 %in% node_colnames) return(match(cand1, node_colnames))
      cand2 <- as.character(k)
      if (cand2 %in% node_colnames) return(match(cand2, node_colnames))
      NA_integer_
    }

    Phi_group_list <- vector("list", n_groups)
    Phi_sets_phi   <- vector("list", n_groups)

    for (g in seq_len(n_groups)) {
      ids <- grp$group_nodes[[g]]
      if (!length(ids)) {
        Phi_group_list[[g]] <- stats::setNames(rep(0, length(spp_all)), spp_all)
        Phi_sets_phi[[g]]   <- character(0)
        next
      }
      col_idx <- vapply(ids, .node_to_col_idx, integer(1L))
      if (anyNA(col_idx)) {
        missing_ids <- ids[is.na(col_idx)]
        stop("`Species.cluster.phi` is missing columns for nodes: ",
             paste(missing_ids, collapse = ", "))
      }
      Phi_g <- Phi[, col_idx, drop = FALSE]
      phi_s <- apply(Phi_g, 1L, max)
      phi_s[phi_s < 0] <- 0

      Phi_group_list[[g]] <- phi_s

      if (strategy %in% c("phi_cover","phi")) {
        Phi_sets_phi[[g]] <- names(phi_s)[phi_s >= min_phi]
      } else {
        Phi_sets_phi[[g]] <- NULL
      }
    }
  }

  ## ---- helper for tie-breaking --------------------------------------------
  .choose_winner <- function(scores_primary,
                             scores_second = NULL,
                             scores_third  = NULL) {
    ok <- !is.na(scores_primary)
    if (!any(ok)) return(NA_character_)

    x1 <- scores_primary[ok]
    max1 <- max(x1)
    if (max1 <= 0) return(NA_character_)

    idx1 <- names(x1)[x1 == max1]
    if (length(idx1) == 1L) return(idx1)

    if (!is.null(scores_second)) {
      s2 <- scores_second[idx1]
      s2_ok <- !is.na(s2)
      if (any(s2_ok)) {
        s2 <- s2[s2_ok]
        max2 <- max(s2)
        idx2 <- names(s2)[s2 == max2]
        if (length(idx2) == 1L) return(idx2)
        idx1 <- idx2
      }
    }

    if (!is.null(scores_third)) {
      s3 <- scores_third[idx1]
      s3_ok <- !is.na(s3)
      if (any(s3_ok)) {
        s3 <- s3[s3_ok]
        max3 <- max(s3)
        idx3 <- names(s3)[s3 == max3]
        if (length(idx3) == 1L) return(idx3)
      }
    }

    "+"
  }

  ## ---- loop over plots ----------------------------------------------------
  for (i in seq_len(N)) {
    pres_sp <- present_list[[i]]

    primary   <- rep(NA_real_, n_groups)
    secondary <- rep(NA_real_, n_groups)
    tertiary  <- rep(NA_real_, n_groups)
    names(primary)   <- group_labels
    names(secondary) <- group_labels
    names(tertiary)  <- group_labels

    for (g in seq_len(n_groups)) {
      if (Ggrp[i, g] <= 0) next

      topo_sp <- intersect(species_sets_topo[[g]], pres_sp)

      if (strategy == "count") {
        n_sp <- length(topo_sp)
        if (n_sp > 0) {
          primary[g] <- n_sp
        }
        next
      }

      if (strategy == "cover") {
        if (!length(topo_sp)) next
        cov_vals <- vm[i, topo_sp, drop = TRUE]
        Cg <- sum(cov_vals)
        if (Cg > 0) {
          primary[g]   <- Cg
          secondary[g] <- length(topo_sp)
        }
        next
      }

      if (strategy == "phi_node") {
        phi_s <- Phi_group_list[[g]]
        if (is.null(phi_s) || !length(phi_s)) next
        if (!length(topo_sp)) next
        phi_use <- phi_s[topo_sp]
        phi_use[is.na(phi_use)] <- 0
        phi_sum <- sum(phi_use)
        if (phi_sum <= 0) next

        cov_vals <- vm[i, topo_sp, drop = TRUE]
        Cg <- sum(cov_vals)
        n_sp <- length(topo_sp)

        primary[g]   <- phi_sum
        secondary[g] <- Cg
        tertiary[g]  <- n_sp
        next
      }

      if (strategy == "phi_cover_node") {
        phi_s <- Phi_group_list[[g]]
        if (is.null(phi_s) || !length(phi_s)) next
        if (!length(topo_sp)) next
        phi_use <- phi_s[topo_sp]
        phi_use[is.na(phi_use)] <- 0
        cov_vals <- vm[i, topo_sp, drop = TRUE]
        score <- sum(cov_vals * phi_use)
        if (score <= 0) next

        n_sp <- length(topo_sp)

        primary[g]   <- score
        secondary[g] <- n_sp
        next
      }

      if (strategy == "phi_cover") {
        phi_s <- Phi_group_list[[g]]
        if (is.null(phi_s) || !length(phi_s)) next
        sp_phi <- Phi_sets_phi[[g]]
        if (!length(sp_phi)) next
        sp_use <- intersect(sp_phi, pres_sp)
        if (!length(sp_use)) next

        phi_use <- phi_s[sp_use]
        phi_use[is.na(phi_use)] <- 0
        cov_vals <- vm[i, sp_use, drop = TRUE]

        score <- sum(cov_vals * phi_use)
        if (score <= 0) next

        Cg <- sum(cov_vals)
        n_sp <- length(sp_use)

        primary[g]   <- score
        secondary[g] <- Cg
        tertiary[g]  <- n_sp
        next
      }

      if (strategy == "phi") {
        phi_s <- Phi_group_list[[g]]
        if (is.null(phi_s) || !length(phi_s)) next
        sp_phi <- Phi_sets_phi[[g]]
        if (!length(sp_phi)) next
        sp_use <- intersect(sp_phi, pres_sp)
        if (!length(sp_use)) next

        phi_use <- phi_s[sp_use]
        phi_use[is.na(phi_use)] <- 0
        phi_sum <- sum(phi_use)
        if (phi_sum <= 0) next

        cov_vals <- vm[i, sp_use, drop = TRUE]
        Cg <- sum(cov_vals)
        n_sp <- length(sp_use)

        primary[g]   <- phi_sum
        secondary[g] <- Cg
        tertiary[g]  <- n_sp
        next
      }
    }

    assigned[i] <- .choose_winner(primary, secondary, tertiary)
  }

  ## ---- collapse small groups ----------------------------------------------
  if (min_group_size > 1L) {
    counts <- table(assigned)
    small  <- names(counts)[counts < min_group_size &
                              !names(counts) %in% c("+", NA_character_)]
    if (length(small)) {
      assigned[assigned %in% small] <- "-"
    }
    collapsed <- small
  } else {
    collapsed <- character(0)
  }

  details <- list(
    strategy         = strategy,
    phi_cut          = phi_cut,
    clusters         = clusters,
    groups_used      = group_labels,
    min_phi          = min_phi,
    min_group_size   = min_group_size,
    collapsed_groups = collapsed
  )

  structure(stats::setNames(assigned, names(assigned)), details = details)
}
