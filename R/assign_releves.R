#' Assign relevés (plots) to Cocktail groups using covers and species-cluster phi
#'
#' @description
#' Assign plots (relevés) to user-defined groups of Cocktail nodes.
#'
#' Groups are defined explicitly by \code{clusters}:
#' \itemize{
#'   \item If \code{clusters} is a vector (e.g. \code{c("c_12","c_27")} or \code{c(12,27)}),
#'         each element defines a separate group.
#'   \item If \code{clusters} is a \strong{list}, each list element defines a \strong{union group}
#'         of several nodes (OR of their plot membership vectors).
#' }
#'
#' Each plot is assigned by scoring all groups and selecting the best-scoring group.
#' Depending on \code{plot_membership}, scoring can be restricted to only those groups
#' where the plot is a Cocktail member (based on \code{x$Plot.cluster > 0}, unioned across
#' nodes within a group), or allowed for all groups.
#'
#' @section Strategies:
#' The \code{strategy} defines how group scores are computed (per plot, per group):
#'
#' \describe{
#'
#' \item{\code{"count"}}{
#' Counts how many of the group's candidate species are present in the plot.
#' The group with the maximum count wins.
#' Ties are labeled \code{"+"}. If all counts are 0, assignment is \code{NA}.
#'
#' Candidate species are always the node-constituting species
#' (from \code{x$Cluster.species}), unioned across nodes within each group.
#' }
#'
#' \item{\code{"cover"}}{
#' Sums the cover values (from \code{vegmatrix}) of the group's candidate species
#' present in the plot. The maximum sum wins.
#' Ties are broken by species count; remaining ties are labeled \code{"+"}.
#' If all sums are 0, assignment is \code{NA}.
#'
#' Candidate species are always the node-constituting species
#' (from \code{x$Cluster.species}), unioned across nodes within each group.
#' }
#'
#' \item{\code{"phi"}}{
#' Requires \code{x$Species.cluster.phi}. For each group and species, a group-level
#' fidelity weight \eqn{\phi(s,g)} is defined as the maximum \eqn{\phi} across the
#' group's nodes; negative values are set to 0.
#'
#' The score is the sum of \eqn{\phi(s,g)} across candidate species present in the plot.
#' Ties are broken by total cover, then by species count; remaining ties are labeled \code{"+"}.
#' If all scores are 0, assignment is \code{NA}.
#'
#' Candidate species depend on \code{min_phi}:
#' \itemize{
#'   \item If \code{min_phi = NULL} (default): candidate species are the node-constituting species
#'         (from \code{x$Cluster.species}).
#'   \item If \code{min_phi} is numeric in \eqn{[0,1]}: candidate species are selected from the
#'         full fidelity profile as those with \eqn{\phi(s,g) \ge min\_phi}.
#' }
#' }
#'
#' \item{\code{"phi_cover"}}{
#' Requires \code{x$Species.cluster.phi}. Uses \eqn{\phi(s,g)} defined as for \code{"phi"},
#' but scores each group as \eqn{\sum \mathrm{cover}(i,s)\,\phi(s,g)} over candidate species
#' present in the plot.
#'
#' Ties are broken by total cover, then by species count; remaining ties are labeled \code{"+"}.
#' If all scores are 0, assignment is \code{NA}.
#'
#' Candidate species depend on \code{min_phi} exactly as in \code{"phi"}.
#' }
#'
#' }
#'
#' @section Missing \code{Species.cluster.phi}:
#' If \code{strategy} is \code{"phi"} or \code{"phi_cover"} but \code{x$Species.cluster.phi}
#' is missing, the function warns and falls back to:
#' \itemize{
#'   \item \code{"count"} if \code{plot_membership = TRUE}, or
#'   \item \code{"cover"} if \code{plot_membership = FALSE}.
#' }
#'
#' @param x A \code{"cocktail"} object (result of \code{\link{cocktail_cluster}}),
#'   containing at least \code{Cluster.species}, \code{Cluster.merged},
#'   \code{Cluster.height}, and \code{Plot.cluster}. For phi-based strategies,
#'   \code{x$Species.cluster.phi} should be present.
#'
#' @param vegmatrix Numeric matrix/data.frame of covers (plots \eqn{\times} species).
#'   Row names must include the plots in \code{x$Plot.cluster}, and column names must
#'   include species used in \code{x$Cluster.species}. Missing values are treated as 0.
#'
#' @param strategy One of \code{c("count","cover","phi","phi_cover")}.
#'
#' @param clusters Mandatory selection of clusters (nodes) defining groups.
#'   Can be:
#'   \itemize{
#'     \item a vector of labels/IDs (e.g. \code{c("c_12","c_27")} or \code{c(12,27)}),
#'           where each element defines a separate group; or
#'     \item a \strong{list} of such vectors, where each element defines a \strong{union group}
#'           of nodes (OR of memberships).
#'   }
#'   Within each group, if both an ancestor and descendant are present,
#'   only the topmost (ancestor) is kept.
#'
#' @param plot_membership Logical. If \code{TRUE} (default), a group competes for a plot
#'   only if the plot is a Cocktail member of that group
#'   (based on \code{x$Plot.cluster > 0}, unioned across nodes in the group).
#'   If \code{FALSE}, all groups compete for all plots.
#'
#' @param min_phi NULL (default) or numeric in \eqn{[0,1]}.
#' Only affects phi-based strategies (\code{"phi"}, \code{"phi_cover"}):
#' \itemize{
#'   \item if \code{NULL}: use node-constituting species only (from \code{x$Cluster.species});
#'   \item if numeric: define candidate species by the full fidelity profile
#'     \eqn{\phi(s,g) \ge min\_phi}.
#' }
#' If \code{min_phi} is not \code{NULL} and is outside \eqn{[0,1]}, an error is raised.
#'
#' @param min_group_size Integer >= 1. After assignment, any group with fewer plots
#'   than this threshold is relabeled \code{"-"}. Default 1.
#'
#' @return
#' A named character vector (names = plot IDs) with values:
#' \itemize{
#'   \item group labels \code{"g_<id>"} or \code{"g_<id>_<id>_..."},
#'   \item \code{"+"} for ties,
#'   \item \code{"-"} for groups collapsed by \code{min_group_size},
#'   \item \code{NA} when no group wins (all scores 0 or no eligible group).
#' }
#' An attribute \code{"details"} is attached with diagnostics.
#'
#' @import Matrix
#' @export
assign_releves <- function(
    x,
    vegmatrix,
    strategy        = c("count","cover","phi","phi_cover"),
    clusters,
    plot_membership = TRUE,
    min_phi         = NULL,
    min_group_size  = 1L
) {

  ## ---- helpers ------------------------------------------------------------
  .is_cocktail <- function(obj) {
    is.list(obj) &&
      all(c("Cluster.species","Cluster.merged","Cluster.height","Plot.cluster") %in% names(obj))
  }

  .parse_clusters_arg <- function(v) {
    if (missing(v) || is.null(v)) return(NULL)

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
      # vector -> each element is its own group
      lapply(as.list(v), parse_one)
    }
  }

  .keep_topmost_within <- function(ids, CM) {
    ids <- sort(unique(ids))
    if (!length(ids)) return(ids)
    kids <- unique(as.integer(CM[ids, , drop = FALSE][CM[ids, , drop = FALSE] > 0]))
    sort(setdiff(ids, intersect(ids, kids)))
  }

  .ensure_pc_colnames <- function(PC) {
    if (is.null(colnames(PC))) {
      colnames(PC) <- paste0("c_", seq_len(ncol(PC)))
    } else {
      num_only <- grepl("^\\d+$", colnames(PC))
      colnames(PC)[num_only] <- paste0("c_", colnames(PC)[num_only])
    }
    PC
  }

  .build_groups <- function(CS, CM, PC, clusters) {
    spp_all <- colnames(CS)
    node_groups <- .parse_clusters_arg(clusters)

    if (is.null(node_groups) || !length(node_groups)) {
      stop("`clusters` must be provided and must define at least one group.")
    }

    # validate IDs + keep only topmost nodes per group
    node_groups <- lapply(node_groups, function(ids) {
      ids <- ids[is.finite(ids)]
      ids <- ids[ids > 0]
      .keep_topmost_within(ids, CM)
    })

    keep_group <- vapply(node_groups, function(v) length(v) > 0, logical(1))
    if (!any(keep_group)) {
      stop("After ancestor/descendant reduction, no groups remain.")
    }
    node_groups <- node_groups[keep_group]

    PC <- .ensure_pc_colnames(PC)

    Glist        <- vector("list", length(node_groups))
    lablist      <- character(length(node_groups))
    species_sets <- vector("list", length(node_groups))

    for (g in seq_along(node_groups)) {
      ids  <- node_groups[[g]]
      cols <- paste0("c_", ids)

      missing_cols <- setdiff(cols, colnames(PC))
      if (length(missing_cols)) {
        stop("Internal mismatch: x$Plot.cluster is missing columns: ",
             paste(missing_cols, collapse = ", "))
      }

      # plot membership: union across nodes in the group
      sub_pc <- PC[, cols, drop = FALSE]
      Gi <- Matrix::rowSums(sub_pc != 0) > 0
      Glist[[g]] <- as.numeric(Gi)

      # group label
      lablist[g] <- if (length(ids) == 1L) {
        paste0("g_", ids)
      } else {
        paste0("g_", paste(ids, collapse = "_"))
      }

      # "membership species set" = union of node-constituting species across nodes
      sp_union_logical <- apply(CS[ids, , drop = FALSE] > 0, 2, any)
      species_sets[[g]] <- spp_all[sp_union_logical]
    }

    G <- do.call(cbind, Glist)
    colnames(G) <- lablist

    list(
      group_nodes            = node_groups,
      group_labels           = lablist,
      G                      = G,
      species_sets_membership = species_sets
    )
  }

  .choose_winner <- function(scores_primary,
                             scores_second = NULL,
                             scores_third  = NULL) {

    ok <- !is.na(scores_primary)
    if (!any(ok)) return(NA_character_)

    x1 <- scores_primary[ok]
    max1 <- max(x1)
    if (!is.finite(max1) || max1 <= 0) return(NA_character_)

    idx1 <- names(x1)[x1 == max1]
    if (length(idx1) == 1L) return(idx1)

    if (!is.null(scores_second)) {
      s2 <- scores_second[idx1]
      s2 <- s2[!is.na(s2)]
      if (length(s2)) {
        max2 <- max(s2)
        idx2 <- names(s2)[s2 == max2]
        if (length(idx2) == 1L) return(idx2)
        idx1 <- idx2
      }
    }

    if (!is.null(scores_third)) {
      s3 <- scores_third[idx1]
      s3 <- s3[!is.na(s3)]
      if (length(s3)) {
        max3 <- max(s3)
        idx3 <- names(s3)[s3 == max3]
        if (length(idx3) == 1L) return(idx3)
      }
    }

    "+"
  }

  ## ---- checks / setup -----------------------------------------------------
  strategy <- match.arg(strategy)

  if (!.is_cocktail(x)) {
    stop("`x` must be a Cocktail object with Cluster.species, Cluster.merged, Cluster.height, and Plot.cluster.")
  }

  if (missing(clusters) || is.null(clusters)) {
    stop("`clusters` must be provided. Use `clusters_at_cut()` to obtain nodes at a cut if needed.")
  }

  if (!is.logical(plot_membership) || length(plot_membership) != 1L || is.na(plot_membership)) {
    stop("`plot_membership` must be a single logical value (TRUE/FALSE).")
  }

  # min_phi validity: NULL or numeric in [0,1]
  if (!is.null(min_phi)) {
    if (!is.numeric(min_phi) || length(min_phi) != 1L || is.na(min_phi)) {
      stop("`min_phi` must be NULL or a single numeric value in [0,1].")
    }
    if (min_phi < 0 || min_phi > 1) {
      stop("`min_phi` must be NULL or a numeric value in [0,1].")
    }
  }

  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  PC <- x$Plot.cluster

  if (!is.matrix(CS)) stop("`x$Cluster.species` must be a matrix.")

  spp_all <- colnames(CS)
  if (is.null(spp_all)) stop("`x$Cluster.species` must have species column names.")

  # vegmatrix
  if (missing(vegmatrix) || is.null(vegmatrix)) {
    stop("`vegmatrix` must be provided (covers).")
  }
  vm <- as.matrix(vegmatrix)
  vm[is.na(vm)] <- 0
  storage.mode(vm) <- "double"

  if (is.null(rownames(vm))) stop("`vegmatrix` must have plot row names.")
  if (is.null(colnames(vm))) stop("`vegmatrix` must have species column names.")

  # plot names from Plot.cluster (or align to vegmatrix if missing)
  plot_names <- rownames(PC)
  if (is.null(plot_names)) {
    plot_names <- rownames(vm)
    rownames(PC) <- plot_names
  }

  # align plots
  if (!all(plot_names %in% rownames(vm))) {
    missing_pl <- setdiff(plot_names, rownames(vm))
    stop("`vegmatrix` is missing plots used in Cocktail: ",
         paste(head(missing_pl, 10), collapse = ", "),
         if (length(missing_pl) > 10) " ..." else "")
  }
  vm <- vm[plot_names, , drop = FALSE]

  # align species
  common <- intersect(colnames(vm), spp_all)
  if (!length(common)) stop("No overlapping species between `vegmatrix` and Cocktail species.")
  vm <- vm[, common, drop = FALSE]
  CS <- CS[, common, drop = FALSE]
  spp_all <- common

  N <- nrow(vm)

  # Ensure Plot.cluster is a Matrix for sparse ops (only used for membership OR)
  if (!inherits(PC, "Matrix")) {
    PC <- Matrix::Matrix(as.matrix(PC), sparse = TRUE)
  }
  PC <- .ensure_pc_colnames(PC)

  # Build group definitions (membership vectors + membership species sets)
  grp <- .build_groups(CS, CM, PC, clusters)
  group_labels            <- grp$group_labels
  Ggrp                    <- as.matrix(grp$G)                   # N x G, 0/1
  species_sets_membership <- grp$species_sets_membership        # list length G
  n_groups                <- ncol(Ggrp)

  min_group_size <- as.integer(min_group_size)
  if (!is.finite(min_group_size) || min_group_size < 1L) min_group_size <- 1L

  assigned <- rep(NA_character_, N)
  names(assigned) <- plot_names

  # present species list per plot
  species_names <- colnames(vm)
  present_list <- lapply(seq_len(N), function(i) species_names[vm[i, ] > 0])

  ## ---- phi precomputation (only when needed) ------------------------------
  need_phi <- strategy %in% c("phi","phi_cover")

  Phi_group_list <- NULL
  Phi_sets_phi   <- NULL

  if (need_phi) {
    if (!("Species.cluster.phi" %in% names(x)) || is.null(x$Species.cluster.phi)) {
      fallback <- if (isTRUE(plot_membership)) "count" else "cover"
      warning(
        "Strategy '", strategy, "' requires x$Species.cluster.phi, which is not present.\n",
        "Recompute cocktail_cluster(..., species_cluster_phi = TRUE) to enable phi-based assignment.\n",
        "Falling back to strategy '", fallback, "'."
      )
      strategy <- fallback
      need_phi <- FALSE
    }
  }

  if (need_phi) {
    Phi <- x$Species.cluster.phi
    if (!is.matrix(Phi)) stop("`x$Species.cluster.phi` must be a matrix.")
    if (is.null(rownames(Phi))) stop("`x$Species.cluster.phi` must have species row names.")
    if (is.null(colnames(Phi))) stop("`x$Species.cluster.phi` must have node column names (e.g. 'c_1').")

    # restrict Phi to overlapping species
    if (!all(spp_all %in% rownames(Phi))) {
      missing_sp <- setdiff(spp_all, rownames(Phi))
      stop("`Species.cluster.phi` is missing species: ",
           paste(head(missing_sp, 10), collapse = ", "),
           if (length(missing_sp) > 10) " ..." else "")
    }
    Phi <- Phi[spp_all, , drop = FALSE]

    node_colnames <- colnames(Phi)

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
      col_idx <- vapply(ids, .node_to_col_idx, integer(1L))

      if (anyNA(col_idx)) {
        stop("`Species.cluster.phi` is missing columns for nodes: ",
             paste(ids[is.na(col_idx)], collapse = ", "))
      }

      # group phi profile: max across nodes
      Phi_g <- Phi[, col_idx, drop = FALSE]
      phi_s <- apply(Phi_g, 1L, max)
      phi_s[phi_s < 0] <- 0
      Phi_group_list[[g]] <- phi_s

      # candidate set for full-profile mode (min_phi numeric)
      if (!is.null(min_phi)) {
        Phi_sets_phi[[g]] <- names(phi_s)[phi_s >= min_phi]
      } else {
        Phi_sets_phi[[g]] <- NULL
      }
    }
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

      # Restrict by Cocktail plot membership if requested
      if (isTRUE(plot_membership) && Ggrp[i, g] <= 0) next

      # Candidate species:
      # - count/cover always use membership species sets (node-constituting)
      # - phi/phi_cover:
      #     * if min_phi NULL -> membership species sets
      #     * if min_phi numeric -> full-profile phi-selected species
      if (strategy %in% c("count","cover") || is.null(min_phi) || !need_phi) {
        sp_use <- intersect(species_sets_membership[[g]], pres_sp)
      } else {
        sp_phi <- Phi_sets_phi[[g]]
        sp_use <- intersect(sp_phi, pres_sp)
      }

      if (!length(sp_use)) next

      if (strategy == "count") {
        primary[g] <- length(sp_use)
        next
      }

      if (strategy == "cover") {
        cov_vals <- vm[i, sp_use, drop = TRUE]
        sc <- sum(cov_vals)
        if (sc > 0) {
          primary[g]   <- sc
          secondary[g] <- length(sp_use)
        }
        next
      }

      if (strategy == "phi") {
        phi_s <- Phi_group_list[[g]]
        phi_use <- phi_s[sp_use]
        phi_use[is.na(phi_use)] <- 0
        sc <- sum(phi_use)
        if (sc > 0) {
          cov_vals <- vm[i, sp_use, drop = TRUE]
          primary[g]   <- sc
          secondary[g] <- sum(cov_vals)
          tertiary[g]  <- length(sp_use)
        }
        next
      }

      if (strategy == "phi_cover") {
        phi_s <- Phi_group_list[[g]]
        phi_use <- phi_s[sp_use]
        phi_use[is.na(phi_use)] <- 0
        cov_vals <- vm[i, sp_use, drop = TRUE]
        sc <- sum(cov_vals * phi_use)
        if (sc > 0) {
          primary[g]   <- sc
          secondary[g] <- sum(cov_vals)
          tertiary[g]  <- length(sp_use)
        }
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
    strategy          = strategy,
    clusters          = clusters,
    plot_membership   = plot_membership,
    min_phi           = min_phi,
    groups_used       = group_labels,
    min_group_size    = min_group_size,
    collapsed_groups  = collapsed
  )

  structure(stats::setNames(assigned, names(assigned)), details = details)
}
