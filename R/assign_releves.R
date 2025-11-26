#' Assign relevés (plots) to Cocktail groups (topology, covers, fuzzy)
#'
#' @description
#' Assign plots to groups derived from a Cocktail tree in four ways:
#'
#' - **"topology"** — Use *topological* membership from `x$Plot.cluster` at a
#'   given `phi_cut` (or user-specified `clusters`). A plot can be a member of
#'   **0, 1, or several** groups. If **every** plot belongs to **exactly one**
#'   group, the function returns a **character vector** of assignments. If **any**
#'   plot has **multiple** memberships, the function returns a **binary matrix**
#'   (plots × groups) with 0/1 membership flags.
#'
#' - **"topology_cover"** — Like *topology*, but when `vegmatrix` contains
#'   cover values, score each plot–group pair by the **sum of transformed covers**
#'   across the group’s *topological* species (i.e., species that define that
#'   parent node). Plots are **eligible** only if they are topological members of
#'   the group (`x$Plot.cluster`) **and** their **raw (pre-transform)** cover sum
#'   for the group's species is at least `min_cover`. The **unique maximum** score
#'   (if any) wins; otherwise `"not assigned"`.
#'
#' - **"topology_fuzzy"** — Compute **direct** (Option A) \eqn{\phi} between each
#'   species and **each selected group** (the group’s plot-membership is the
#'   OR-union of its parent nodes’ memberships). Then **mask** these \eqn{\phi}
#'   weights so that only the group’s *topological species* (from the Cocktail
#'   node) can contribute. Plot scores are (transformed covers) × (masked \eqn{\phi}).
#'   Eligibility by `min_cover` uses the group’s *topological* species. Unique
#'   maximum rule as above.
#'
#' - **"fuzzy"** — Compute **direct** (Option A) \eqn{\phi} for species against
#'   **each selected group**, then **filter species by fidelity**: keep species
#'   with \eqn{\phi \ge} `phi_fuzzy` for that group (`phi_fuzzy` defaults to
#'   `phi_cut`). Plot scores are (transformed covers) × (filtered \eqn{\phi}).
#'   Eligibility by `min_cover` uses the **filtered species set** per group.
#'   Unique maximum rule as above.
#'
#' @param x A Cocktail object (list) with components
#'   `Cluster.species`, `Cluster.merged`, `Cluster.height`, and `Plot.cluster`
#'   (e.g., from `cocktail_cluster()`).
#' @param vegmatrix Optional numeric matrix/data.frame of **covers** with plots
#'   in rows and species in columns. Needed for `"topology_cover"`,
#'   `"topology_fuzzy"`, and `"fuzzy"`. Missing values are treated as 0.
#' @param strategy One of `c("topology","topology_cover","topology_fuzzy","fuzzy")`.
#' @param phi_cut Numeric in (0,1). Height cut used to choose **parent clusters**
#'   (unless `clusters` is supplied).
#' @param clusters Optional selection of clusters **instead of** `phi_cut`.
#'   - If a **vector** of labels/IDs (e.g., `c("c_123","c_456")` or `c(123,456)`),
#'     each element forms a separate group.
#'   - If a **list**, each element is a **union** (group) of several nodes,
#'     e.g., `list(c("c_123","c_234"), 567)`.
#'   Within each group, if both an **ancestor** and its **descendant** are
#'   present, the function **keeps only the topmost** (ancestor) to avoid
#'   double-counting.
#' @param cover_transform Cover transform: `c("sqrt","none","log1p")`. Default `"sqrt"`.
#' @param min_cover Numeric percent in (0,100). A plot–group pair is **ineligible**
#'   unless the **raw (pre-transform)** sum of the group’s species in the plot is
#'   at least this value. Default 0.
#' @param phi_fuzzy Fidelity threshold for `"fuzzy"` (defaults to `phi_cut`).
#'   Species with \eqn{\phi <} `phi_fuzzy` are dropped for that group.
#' @param min_group_size Integer ≥ 1. After assignment, any group with fewer than
#'   this many plots is relabeled `"not assigned"`. Default 1 (no collapsing).
#'
#' @details
#' **Parent clusters** at a cut are the “topmost” nodes with height ≥ `phi_cut`
#' that **do not** have an ancestor also at/above the cut. If `clusters` is not
#' provided, each parent cluster at the cut forms **one group**. If `clusters`
#' is a **list**, each element defines a **union** (OR) of several nodes; plot
#' membership for that group is the OR across those nodes’ memberships in
#' `x$Plot.cluster`, and **species sets for topology** are the **unions** of the
#' nodes’ topological species (`x$Cluster.species` at those node IDs).
#'
#' **Fuzzy (Option A, direct):** for species–group association we *recompute*
#' \eqn{\phi} directly against the **group** membership vector (the OR-union
#' of its components), i.e., using the 2×2 counts \eqn{a,b,c,d} per species vs.
#' the group’s plot membership. This yields diagnostics for the **union as a
#' single unit** (statistically exact for the redefined group).
#'
#' **Topology output shape:** In `"topology"`, a plot can belong to multiple
#' groups. If **all** plots have exactly one membership, a **vector** of labels
#' is returned; otherwise a **binary matrix** (plots × groups) is returned.
#'
#' @return
#' - For `"topology"`: character vector of assignments **or** a binary membership
#'   matrix (plots × groups) if multiple memberships occur.
#' - For the other strategies: a character vector of assigned group labels
#'   (or `"not assigned"`).
#' The object carries an attribute `"details"` with inputs and diagnostics
#' (e.g., scores matrix, groups used, species sets, etc.).
#'
#' @seealso \code{\link{clusters_at_cut}}
#' @export
assign_releves <- function(
    x,
    vegmatrix       = NULL,
    strategy        = c("topology","topology_cover","topology_fuzzy","fuzzy"),
    phi_cut         = 0.30,
    clusters        = NULL,     # vector or list (for unions)
    cover_transform = c("sqrt","none","log1p"),
    min_cover       = 0,
    phi_fuzzy       = NULL,     # defaults to phi_cut
    min_group_size  = 1
) {
  ## ---- helpers ------------------------------------------------------------
  .is_cocktail <- function(obj) {
    is.list(obj) && all(c("Cluster.species","Cluster.merged","Cluster.height","Plot.cluster") %in% names(obj))
  }
  .cover_transform <- function(xx, how) {
    if (how == "sqrt") sqrt(pmax(xx, 0))
    else if (how == "log1p") log1p(pmax(xx, 0))
    else xx
  }
  # parents (topmost) at a cut
  .parents_at_cut <- function(CM, H, cut) {
    idx <- which(H >= cut); if (!length(idx)) return(integer(0))
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
  # build groups from phi_cut or clusters arg: returns list with slots
  # $group_nodes (list of integer vectors), $group_labels, $G (N×G membership), $species_sets (list)
  .build_groups <- function(CS, CM, H, PC, phi_cut, clusters) {
    spp_all <- colnames(CS)
    N <- nrow(PC)
    # choose node sets per group
    node_groups <- if (is.null(clusters)) {
      as.list(.parents_at_cut(CM, H, phi_cut))
    } else {
      .parse_clusters_arg(clusters)
    }
    if (!length(node_groups)) stop("No groups available (check phi_cut or clusters).")
    # validate & reduce each group to topmost
    node_groups <- lapply(seq_along(node_groups), function(i) {
      ids <- node_groups[[i]]
      if (any(!is.finite(ids))) stop("`clusters` contains non-integer/invalid IDs.")
      ids <- ids[ids > 0]
      .keep_topmost_within(ids, CM)
    })
    # discard empty groups (after reduction)
    keep_group <- vapply(node_groups, function(v) length(v) > 0, logical(1))
    if (!any(keep_group)) stop("After ancestor/descendant reduction, no groups remain.")
    node_groups <- node_groups[keep_group]

    # membership per group: OR of PC[, paste0("c_", id)]
    Glist <- list()
    lablist <- character(length(node_groups))
    species_sets <- vector("list", length(node_groups))
    for (g in seq_along(node_groups)) {
      ids <- node_groups[[g]]
      cols <- paste0("c_", ids)
      missing_cols <- setdiff(cols, colnames(PC))
      if (length(missing_cols))
        stop("Internal mismatch: Plot.cluster is missing: ", paste(missing_cols, collapse = ", "))
      Gi <- rowSums(PC[, cols, drop = FALSE] > 0) > 0
      Glist[[g]] <- as.numeric(Gi)              # 0/1 numeric
      lablist[g] <- paste0("g_", paste(ids, collapse = "_"))
      # union of topological species from these node IDs
      sp_union <- colnames(CS)[apply(CS[ids, , drop = FALSE] > 0, 2, any)]
      species_sets[[g]] <- sp_union
    }
    G <- do.call(cbind, Glist); colnames(G) <- lablist; rownames(G) <- rownames(PC)
    list(group_nodes = node_groups, group_labels = lablist, G = G, species_sets = species_sets)
  }
  # Direct (Option A) phi for species vs groups given X (N×S pres/abs) and G (N×G group membership)
  .phi_direct_species_by_groups <- function(X, G) {
    # X: N×S (logical/numeric 0/1), G: N×G (0/1)
    storage.mode(X) <- "double"; storage.mode(G) <- "double"
    N <- nrow(X); S <- ncol(X); Gk <- ncol(G)
    a  <- crossprod(X, G)                  # S×G
    p  <- colSums(X)                       # length S
    g1 <- colSums(G)                       # length G
    b <- matrix(p,  nrow = S, ncol = Gk) - a
    c <- matrix(g1, nrow = S, ncol = Gk, byrow = TRUE) - a
    d <- (N - matrix(g1, nrow = S, ncol = Gk, byrow = TRUE)) - b
    den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
    Phi <- (a * d - b * c) / den
    Phi[!is.finite(den) | den <= 0] <- 0
    Phi
  }

  ## ---- checks / setup -----------------------------------------------------
  strategy        <- match.arg(strategy)
  cover_transform <- match.arg(cover_transform)

  if (!.is_cocktail(x)) {
    stop("x must be a Cocktail object with Cluster.species, Cluster.merged, Cluster.height, and Plot.cluster.")
  }
  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  PC <- x$Plot.cluster

  ## ---- normalize Plot.cluster colnames to "c_<index>" ---------
  if (is.null(colnames(PC))) {
    colnames(PC) <- paste0("c_", seq_len(ncol(PC)))
  } else {
    num_only <- grepl("^\\d+$", colnames(PC))
    colnames(PC)[num_only] <- paste0("c_", colnames(PC)[num_only])
  }

  plot_names <- rownames(PC); if (is.null(plot_names)) plot_names <- paste0("plot_", seq_len(nrow(PC)))
  spp_all <- colnames(CS)

  # build groups (nodes, labels, OR-membership, topology species sets)
  grp <- .build_groups(CS, CM, H, PC, phi_cut, clusters)
  group_labels <- grp$group_labels
  Ggrp <- as.matrix(grp$G)                   # N × G
  species_sets_topo <- grp$species_sets      # length G (character vectors)

  # covers if needed
  need_cover <- strategy %in% c("topology_cover","topology_fuzzy","fuzzy")
  if (need_cover) {
    if (is.null(vegmatrix)) stop("`vegmatrix` is required for strategy = '", strategy, "'.")
    vm <- as.matrix(vegmatrix); vm[is.na(vm)] <- 0
    storage.mode(vm) <- "double"
    if (is.null(colnames(vm))) stop("`vegmatrix` must have species column names.")
    if (is.null(rownames(vm))) rownames(vm) <- plot_names
    # align species
    common <- intersect(colnames(vm), spp_all)
    if (!length(common)) stop("No overlapping species between covers and Cocktail species.")
    vm <- vm[, common, drop = FALSE]
    # light diagnostics
    if (all(vm %in% c(0,1))) {
      message("Covers appear binary; results will be presence-like.")
    } else if (length(unique(as.numeric(vm))) <= 10) {
      message("Covers look ordinal/limited; percentage covers are recommended for best performance.")
    }
  }

  # hygiene
  min_cover <- max(0, min(100, as.numeric(min_cover)))
  min_group_size <- as.integer(min_group_size)
  if (!is.finite(min_group_size) || min_group_size < 1) min_group_size <- 1L
  if (is.null(phi_fuzzy)) phi_fuzzy <- phi_cut
  phi_fuzzy <- max(0, min(1, as.numeric(phi_fuzzy)))

  # output scaffolds
  assigned <- rep("not assigned", nrow(Ggrp)); names(assigned) <- plot_names
  details <- list(
    strategy        = strategy,
    phi_cut         = phi_cut,
    phi_fuzzy       = phi_fuzzy,
    groups_used     = group_labels,
    min_cover       = min_cover,
    cover_transform = cover_transform,
    min_group_size  = min_group_size
  )

  ## ---- strategies ---------------------------------------------------------
  if (strategy == "topology") {
    hits <- (Ggrp > 0)           # N × G
    per_row <- rowSums(hits)
    if (all(per_row %in% c(0,1))) {
      # unique or none → vector
      winner <- max.col(hits, ties.method = "first")
      ok <- (per_row == 1)
      assigned[ok] <- group_labels[winner[ok]]
      # collapse tiny groups if requested
      if (min_group_size > 1) {
        counts <- table(assigned)
        small  <- names(counts)[counts < min_group_size & names(counts) != "not assigned"]
        if (length(small)) assigned[assigned %in% small] <- "not assigned"
        details$collapsed_groups <- small
      } else details$collapsed_groups <- character(0)
      structure(stats::setNames(assigned, names(assigned)), details = details)
    } else {
      # multi-membership exists → return binary matrix
      M <- hits * 1L
      rownames(M) <- plot_names; colnames(M) <- group_labels
      attr(M, "details") <- details
      return(M)
    }

  } else if (strategy == "topology_cover") {
    # species sets = topology sets
    W <- .cover_transform(vm, cover_transform)      # transformed covers
    # raw sums per group for eligibility
    S_raw <- matrix(0, nrow = nrow(vm), ncol = ncol(Ggrp),
                    dimnames = list(rownames(vm), group_labels))
    S <- S_raw
    for (g in seq_along(species_sets_topo)) {
      sp <- intersect(species_sets_topo[[g]], colnames(vm))
      if (length(sp)) {
        S_raw[, g] <- rowSums(vm[, sp, drop = FALSE])
        S[, g]     <- rowSums(W[,  sp, drop = FALSE])
      }
    }
    eligible <- (Ggrp > 0) & (S_raw >= min_cover)
    S[!eligible] <- -Inf

    winner <- max.col(S, ties.method = "first")
    best   <- S[cbind(seq_len(nrow(S)), winner)]
    ties   <- apply(S, 1, function(v) sum(v == max(v)))
    ok <- is.finite(best) & (ties == 1)
    assigned[ok] <- group_labels[winner[ok]]

    if (min_group_size > 1) {
      counts <- table(assigned)
      small  <- names(counts)[counts < min_group_size & names(counts) != "not assigned"]
      if (length(small)) assigned[assigned %in% small] <- "not assigned"
      details$collapsed_groups <- small
    } else details$collapsed_groups <- character(0)

    details$scores <- S
    return(structure(stats::setNames(assigned, names(assigned)), details = details))

  } else {
    # --- fuzzy family: direct (Option A) phi vs groups ---
    # presence/absence for species
    X <- (vm > 0); storage.mode(X) <- "double"   # N × S'
    # recompute phi S' × G directly
    Phi <- .phi_direct_species_by_groups(X, Ggrp)
    colnames(Phi) <- group_labels
    rownames(Phi) <- colnames(X)

    if (strategy == "topology_fuzzy") {
      # mask Phi to only topology species per group
      for (g in seq_along(species_sets_topo)) {
        keep_sp <- intersect(species_sets_topo[[g]], rownames(Phi))
        mask <- !(rownames(Phi) %in% keep_sp)
        Phi[mask, g] <- 0
      }
      # eligibility set = topology species (raw covers)
      S_elig_raw <- matrix(0, nrow = nrow(vm), ncol = ncol(Ggrp),
                           dimnames = list(rownames(vm), group_labels))
      for (g in seq_along(species_sets_topo)) {
        sp <- intersect(species_sets_topo[[g]], colnames(vm))
        if (length(sp)) S_elig_raw[, g] <- rowSums(vm[, sp, drop = FALSE])
      }
    } else { # strategy == "fuzzy"
      # drop species by phi_fuzzy threshold per group
      Phi[Phi < phi_fuzzy] <- 0
      # eligibility set = species kept by phi_fuzzy (raw covers)
      S_elig_raw <- matrix(0, nrow = nrow(vm), ncol = ncol(Ggrp),
                           dimnames = list(rownames(vm), group_labels))
      for (g in seq_along(group_labels)) {
        sp_keep <- rownames(Phi)[Phi[, g] > 0]
        if (length(sp_keep)) S_elig_raw[, g] <- rowSums(vm[, sp_keep, drop = FALSE])
      }
    }

    # Scores = transformed covers × Phi (S' × G)
    W <- .cover_transform(vm, cover_transform)
    S <- as.matrix(W %*% Phi)
    # eligibility: min_cover only (no topology requirement for fuzzy)
    S[S_elig_raw < min_cover] <- -Inf

    winner <- max.col(S, ties.method = "first")
    best   <- S[cbind(seq_len(nrow(S)), winner)]
    ties   <- apply(S, 1, function(v) sum(v == max(v)))
    ok <- is.finite(best) & (ties == 1)
    assigned[ok] <- group_labels[winner[ok]]

    if (min_group_size > 1) {
      counts <- table(assigned)
      small  <- names(counts)[counts < min_group_size & names(counts) != "not assigned"]
      if (length(small)) assigned[assigned %in% small] <- "not assigned"
      details$collapsed_groups <- small
    } else details$collapsed_groups <- character(0)

    details$scores    <- S
    details$phi_matrix <- Phi
    return(structure(stats::setNames(assigned, names(assigned)), details = details))
  }
}
