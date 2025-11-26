#' Hierarchical plot assignment across multiple phi levels
#'
#' @description
#' Assign each plot to a Cocktail cluster **at all φ levels >= `phi_cut`** and
#' keep a single final assignment (by default from the **finest** level, i.e.,
#' the highest φ where the plot has a unique winner). Supports three strategies:
#'
#' - **"cover"** — like "topology_cover" (scores by transformed covers of a
#'   group's *topological* species), but computed at each φ level.
#' - **"topology_fuzzy"** — like in \code{assign_releves()}:
#'   recompute **direct (Option A)** φ for species vs. the group (OR-union of
#'   the group's nodes; here each group is a single parent node), then **mask**
#'   to *topological species* for that node; score by covers × masked φ.
#' - **"fuzzy"** — recompute **direct** φ for species vs. group; then drop
#'   species with φ < `phi_fuzzy`; score by covers × filtered φ.
#'
#' For each level we enforce **unique maximum** per plot among eligible groups.
#' Eligibility requires (raw, pre-transform) cover sum for the group's species
#' ≥ `min_cover`. (For "fuzzy", eligibility uses the **filtered** species set.)
#'
#' @param x Cocktail object (list) with components:
#'   `Cluster.species`, `Cluster.merged`, `Cluster.height`, `Plot.cluster`.
#' @param tab Numeric matrix/data.frame of **covers** (plots × species).
#'   Required for all strategies here.
#' @param strategy One of `c("cover","topology_fuzzy","fuzzy")`.
#' @param phi_cut Numeric in (0,1). Only φ levels **≥ phi_cut** are considered.
#' @param cover_transform Cover transform: `c("sqrt","none","log1p")`. Default `"sqrt"`.
#' @param min_cover Numeric percent in (0, 100). Plot–group pair is ineligible
#'   unless raw (pre-transform) sum of the group's species in the plot ≥ this value.
#' @param phi_fuzzy Fidelity cutoff for `"fuzzy"` (defaults to `phi_cut`).
#' @param choose_level Which level to use for the final single assignment:
#'   `"finest"` (default; the highest φ with unique assignment), or `"coarsest"`
#'   (the lowest φ ≥ `phi_cut` with unique assignment). You can always switch
#'   later using the `levels$level_id` information stored in the result.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{assignment}: named character vector (plots → group label at the chosen level),
#'         or `"not assigned"` when no unique winner at any level.
#'   \item \code{per_level}: data.frame of per-level assignments (long form):
#'         plot, level_id, phi, group (or "not assigned").
#'   \item \code{levels}: data.frame describing all groups at all levels with
#'         \code{level_id}, \code{phi}, \code{group}, \code{node_id},
#'         \code{parent_group} (in the next coarser level), \code{n_plots}.
#'   \item \code{details}: list with inputs and diagnostics (e.g., transforms,
#'         strategy, and a small map to help remap to other levels later).
#' }
#'
#' @export
assign_releves_new <- function(
    x,
    tab,
    strategy        = c("cover","topology_fuzzy","fuzzy"),
    phi_cut         = 0.30,
    cover_transform = c("sqrt","none","log1p"),
    min_cover       = 0,
    phi_fuzzy       = NULL,
    choose_level    = c("finest","coarsest")
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
  # parents (topmost) at a cut (keep nodes with H >= cut that have no ancestor with H >= cut)
  .parents_at_cut <- function(CM, H, cut) {
    idx <- which(H >= cut)
    if (!length(idx)) return(integer(0))
    kids <- CM[idx, , drop = FALSE]
    children <- unique(as.integer(kids[kids > 0]))
    sort(setdiff(idx, intersect(idx, children)))
  }
  # ancestor relation utilities --------------------------------------------
  .parent_of <- function(node, CM) {
    # CM has columns c("left","right") (standard Cocktail merges): children -> parent relation is implicit.
    # Here we need to find the *parent* of 'node': scan rows where node appears as a child.
    hits <- which(CM[,1] == node | CM[,2] == node)
    if (!length(hits)) return(NA_integer_)
    # In Cocktail trees a node can have at most one parent in the binary merge structure.
    hits[1]
  }
  .is_ancestor <- function(anc, desc, CM) {
    if (anc == desc) return(TRUE)
    cur <- desc
    while (is.finite(cur) && !is.na(cur)) {
      p <- .parent_of(cur, CM)
      if (is.na(p)) return(FALSE)
      if (p == anc) return(TRUE)
      cur <- p
    }
    FALSE
  }
  .map_child_to_parent_group <- function(child_nodes, parent_nodes, CM) {
    # For each child node id, pick the (single) parent group node among parent_nodes that is its ancestor.
    sapply(child_nodes, function(nd) {
      hits <- parent_nodes[vapply(parent_nodes, function(p) .is_ancestor(p, nd, CM), logical(1))]
      if (!length(hits)) NA_integer_ else hits[1]
    })
  }
  # φ for species vs groups given binary X (N×S) and G (N×G)
  .phi_direct_species_by_groups <- function(X, G) {
    storage.mode(X) <- "double"; storage.mode(G) <- "double"
    N <- nrow(X); S <- ncol(X); K <- ncol(G)
    a  <- crossprod(X, G)                       # S×K
    p  <- colSums(X)                            # length S
    g1 <- colSums(G)                            # length K
    b <- matrix(p,  nrow = S, ncol = K) - a
    c <- matrix(g1, nrow = S, ncol = K, byrow = TRUE) - a
    d <- (N - matrix(g1, nrow = S, ncol = K, byrow = TRUE)) - b
    den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
    Phi <- (a * d - b * c) / den
    Phi[!is.finite(den) | den <= 0] <- 0
    Phi
  }

  ## ---- checks / setup -----------------------------------------------------
  strategy        <- match.arg(strategy)
  choose_level    <- match.arg(choose_level)
  cover_transform <- match.arg(cover_transform)

  if (!.is_cocktail(x)) {
    stop("x must be a Cocktail object with Cluster.species, Cluster.merged, Cluster.height, and Plot.cluster.")
  }
  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  PC <- x$Plot.cluster

  # --- normalize Plot.cluster colnames to "c_<index>" (IMPORTANT) ---
  if (is.null(colnames(PC))) {
    colnames(PC) <- paste0("c_", seq_len(ncol(PC)))
  } else {
    num_only <- grepl("^\\d+$", colnames(PC))
    colnames(PC)[num_only] <- paste0("c_", colnames(PC)[num_only])
  }

  # covers (required here)
  if (is.null(tab)) stop("`tab` (covers) must be provided.")
  vm <- as.matrix(tab); vm[is.na(vm)] <- 0
  storage.mode(vm) <- "double"
  if (is.null(colnames(vm))) stop("`tab` must have species column names.")
  if (is.null(rownames(vm))) rownames(vm) <- rownames(PC)
  # align species to Cocktail species set
  spp_all <- colnames(CS)
  common  <- intersect(colnames(vm), spp_all)
  if (!length(common)) stop("No overlapping species between covers and Cocktail species.")
  vm <- vm[, common, drop = FALSE]

  if (is.null(phi_fuzzy)) phi_fuzzy <- phi_cut
  phi_fuzzy <- max(0, min(1, as.numeric(phi_fuzzy)))
  min_cover <- max(0, min(100, as.numeric(min_cover)))

  plot_names <- rownames(PC); if (is.null(plot_names)) plot_names <- paste0("plot_", seq_len(nrow(PC)))

  ## ---- build φ-levels (>= phi_cut), descending (finest first) ------------
  phi_levels <- sort(unique(as.numeric(H[is.finite(H) & H >= phi_cut])), decreasing = TRUE)
  if (!length(phi_levels)) stop("No nodes with height >= phi_cut = ", phi_cut)

  # data containers
  per_level_rows <- list()
  level_desc_rows <- list()

  # For mapping parents between adjacent levels:
  prev_level_nodes <- NULL
  prev_level_labels <- NULL
  prev_level_id <- NULL

  level_id_counter <- 0L

  # precompute transformed covers
  W <- .cover_transform(vm, cover_transform)

  ## ---- iterate levels -----------------------------------------------------
  for (phi in phi_levels) {
    level_id_counter <- level_id_counter + 1L
    level_id <- paste0("L", level_id_counter)

    # parent nodes (topmost) at this cut
    nodes <- .parents_at_cut(CM, H, phi)
    if (!length(nodes)) next

    # groups for this level: each node is one group
    group_nodes  <- as.integer(nodes)
    group_labels <- paste0("g_", group_nodes)
    # membership per group from Plot.cluster
    # ensure columns exist
    cols <- paste0("c_", group_nodes)
    missing_cols <- setdiff(cols, colnames(PC))
    if (length(missing_cols)) stop("Internal mismatch: Plot.cluster is missing: ", paste(missing_cols, collapse = ", "))
    G <- as.matrix(PC[, cols, drop = FALSE] > 0) * 1.0
    colnames(G) <- group_labels
    rownames(G) <- plot_names

    # topology species per group (from Cluster.species rows = nodes)
    topo_species <- lapply(group_nodes, function(id) colnames(CS)[CS[id, , drop = FALSE] > 0])
    names(topo_species) <- group_labels

    # ----- scoring & assignment at this level -----
    # Build eligibility matrix based on raw covers of the group's species
    S_elig_raw <- matrix(0, nrow = nrow(vm), ncol = ncol(G),
                         dimnames = list(rownames(vm), group_labels))

    if (strategy == "cover" || strategy == "topology_fuzzy") {
      for (g in seq_along(group_labels)) {
        sp <- intersect(topo_species[[g]], colnames(vm))
        if (length(sp)) S_elig_raw[, g] <- rowSums(vm[, sp, drop = FALSE])
      }
    }

    if (strategy == "cover") {
      # score = sum of transformed covers across topology species
      S <- matrix(0, nrow = nrow(vm), ncol = ncol(G),
                  dimnames = list(rownames(vm), group_labels))
      for (g in seq_along(group_labels)) {
        sp <- intersect(topo_species[[g]], colnames(vm))
        if (length(sp)) S[, g] <- rowSums(W[, sp, drop = FALSE])
      }
      # eligibility ALSO requires topological membership (plot ∈ group) & min_cover
      eligible <- (G > 0) & (S_elig_raw >= min_cover)
      S[!eligible] <- -Inf

    } else {
      # fuzzy family: recompute direct φ against the group membership
      X   <- (vm > 0) * 1.0
      Phi <- .phi_direct_species_by_groups(X, G)     # S' × K
      colnames(Phi) <- group_labels
      rownames(Phi) <- colnames(X)

      if (strategy == "topology_fuzzy") {
        # mask φ to topology species only
        for (g in seq_along(group_labels)) {
          keep <- intersect(topo_species[[g]], rownames(Phi))
          if (length(keep)) {
            mask <- !(rownames(Phi) %in% keep)
            Phi[mask, g] <- 0
          } else {
            # if no species, whole column will be zeros
            Phi[, g] <- 0
          }
        }
        # score = W %*% Phi; eligibility = min_cover based on topology species
        S <- as.matrix(W %*% Phi)
        S[S_elig_raw < min_cover] <- -Inf

      } else { # strategy == "fuzzy"
        Phi[Phi < phi_fuzzy] <- 0
        # eligibility uses the **kept** species set per group
        S_elig_raw[] <- 0
        for (g in seq_along(group_labels)) {
          keep_sp <- rownames(Phi)[Phi[, g] > 0]
          if (length(keep_sp)) S_elig_raw[, g] <- rowSums(vm[, keep_sp, drop = FALSE])
        }
        S <- as.matrix(W %*% Phi)
        S[S_elig_raw < min_cover] <- -Inf
      }
    }

    # unique max per plot
    winner <- max.col(S, ties.method = "first")
    best   <- S[cbind(seq_len(nrow(S)), winner)]
    ties   <- apply(S, 1, function(v) sum(v == max(v)))
    ok     <- is.finite(best) & (ties == 1)

    assign_vec <- rep("not assigned", nrow(S))
    names(assign_vec) <- rownames(S)
    assign_vec[ok] <- colnames(S)[winner[ok]]

    # stash per-level assignments (long form)
    per_level_rows[[length(per_level_rows) + 1L]] <- data.frame(
      plot     = rownames(S),
      level_id = level_id,
      phi      = phi,
      group    = unname(assign_vec),
      stringsAsFactors = FALSE
    )

    # level descriptor table
    # map each current node to a parent group (if previous level existed and is coarser)
    parent_group <- rep(NA_character_, length(group_nodes))
    if (!is.null(prev_level_nodes)) {
      # previous iteration had higher phi (finer). For hierarchy, we want to connect
      # current (finer) groups to the *next coarser* level. Since we iterate from finest→coarser
      # we assign parents from previous (finer) to current (coarser) in the next loop.
      # So do nothing here; we'll fill this when we move to the next (coarser) level below.
    }

    level_desc_rows[[length(level_desc_rows) + 1L]] <- data.frame(
      level_id = level_id,
      phi      = phi,
      group    = group_labels,
      node_id  = group_nodes,
      parent_group = NA_character_,   # filled when we step to the *next* coarser level
      n_plots  = colSums(G),
      stringsAsFactors = FALSE
    )

    # carry forward for parent mapping (we’ll link on the next loop step)
    prev_level_nodes  <- group_nodes
    prev_level_labels <- group_labels
    prev_level_id     <- level_id
  }

  # We iterated finest→coarser. To fill parent relations, we need to traverse
  # again from fine to coarse and connect *fine* groups to the *next coarser* level.
  # We'll reconstruct a combined table and link by ancestry.
  levels_df <- do.call(rbind, level_desc_rows)
  # order by phi descending (finest first)
  levels_df <- levels_df[order(levels_df$phi, decreasing = TRUE), , drop = FALSE]
  levels_df$parent_group <- NA_character_

  # create lookups by level order
  ord_levels <- unique(levels_df$level_id)  # already finest→coarser
  for (i in seq_along(ord_levels)[-length(ord_levels)]) {
    fine_id   <- ord_levels[i]
    coarse_id <- ord_levels[i + 1L]
    fine_rows   <- which(levels_df$level_id == fine_id)
    coarse_rows <- which(levels_df$level_id == coarse_id)
    fine_nodes   <- levels_df$node_id[fine_rows]
    coarse_nodes <- levels_df$node_id[coarse_rows]
    # for each fine node, find its ancestor among coarse_nodes
    parent_ids <- .map_child_to_parent_group(fine_nodes, coarse_nodes, CM)
    # translate node ids to group labels of the coarse level
    parent_labs <- ifelse(is.na(parent_ids), NA_character_,
                          paste0("g_", parent_ids))
    levels_df$parent_group[fine_rows] <- parent_labs
  }

  # Combine per-level assignments
  per_level <- do.call(rbind, per_level_rows)
  per_level$plot <- as.character(per_level$plot)

  # Final single assignment per plot
  # choose "finest": the first (highest φ) level where group != "not assigned"
  # choose "coarsest": the last (lowest φ) level where group != "not assigned"
  if (choose_level == "finest") {
    # per plot, take the first non-"not assigned" by φ desc
    per_level <- per_level[order(per_level$phi, decreasing = TRUE), , drop = FALSE]
    pick <- tapply(seq_len(nrow(per_level)),
                   per_level$plot,
                   function(ix) {
                     j <- ix[which(per_level$group[ix] != "not assigned")[1]]
                     if (length(j)) j else NA_integer_
                   })
  } else {
    # coarsest: take last non-"not assigned" by φ desc (i.e., minimum φ)
    per_level <- per_level[order(per_level$phi, decreasing = TRUE), , drop = FALSE]
    pick <- tapply(seq_len(nrow(per_level)),
                   per_level$plot,
                   function(ix) {
                     ok_ix <- ix[per_level$group[ix] != "not assigned"]
                     j <- if (length(ok_ix)) ok_ix[length(ok_ix)] else NA_integer_
                     if (length(j)) j else NA_integer_
                   })
  }

  final_assign <- rep("not assigned", length(unique(per_level$plot)))
  names(final_assign) <- sort(unique(per_level$plot))
  sel <- as.integer(pick[ names(final_assign) ])
  good <- is.finite(sel) & !is.na(sel)
  final_assign[good] <- per_level$group[sel[good]]

  # Output
  details <- list(
    strategy        = strategy,
    phi_cut         = phi_cut,
    phi_levels      = phi_levels,
    cover_transform = cover_transform,
    min_cover       = min_cover,
    phi_fuzzy       = phi_fuzzy,
    choose_level    = choose_level
  )

  out <- list(
    assignment = final_assign,
    per_level  = per_level[order(per_level$plot, -per_level$phi), , drop = FALSE],
    levels     = levels_df,
    details    = details
  )
  class(out) <- c("assign_releves_hier", class(out))
  out
}
