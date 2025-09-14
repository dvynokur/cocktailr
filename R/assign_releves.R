# --- helpers ---------------------------------------------------------------

.cover_transform <- function(x, how = c("sqrt", "none", "log1p")) {
  how <- match.arg(how)
  if (how == "sqrt") sqrt(pmax(x, 0))
  else if (how == "log1p") log1p(pmax(x, 0))
  else x
}

# Compare strategy used by both functions.
# scores: matrix [plots x groups]
.pick_groups <- function(scores, compare = c("max", "one_vs_rest", "pairwise")) {
  compare <- match.arg(compare)
  rn <- rownames(scores); cn <- colnames(scores)
  if (is.null(rn)) rn <- seq_len(nrow(scores))
  if (is.null(cn)) cn <- paste0("G", seq_len(ncol(scores)))

  best_idx <- max.col(scores, ties.method = "first")
  best <- cn[best_idx]
  best_score <- scores[cbind(seq_len(nrow(scores)), best_idx)]

  second_score <- apply(scores, 1, function(v) {
    if (length(v) == 1) return(0)
    sort(v, decreasing = TRUE)[pmin(2, length(v))]
  })
  margin <- best_score - second_score

  if (compare == "max") {
    assigned <- best
    stat <- margin
  } else if (compare == "one_vs_rest") {
    sum_all <- rowSums(scores)
    rest <- pmax(sum_all - best_score, 0)
    ratio <- ifelse(rest > 0, best_score / rest, Inf)
    assigned <- best
    stat <- ratio
  } else { # pairwise: must beat every other group
    assigned <- best
    for (i in seq_len(nrow(scores))) {
      b <- best_idx[i]
      if (any(scores[i, -b] >= scores[i, b])) assigned[i] <- "tie"
    }
    stat <- margin
  }

  data.frame(
    plot = rn,
    assigned_group = assigned,
    best_score = best_score,
    second_best = second_score,
    margin = margin,
    compare_stat = stat,
    stringsAsFactors = FALSE
  )
}

# --- 1) Strict assignment at a phi cut -------------------------------------

#' Assign releves to strict Cocktail groups at a phi cut
#'
#' @param x Result of `cocktail_cluster()`.
#' @param vegmatrix numeric matrix/data.frame: plots x species (covers).
#' @param phi numeric; cut level (e.g. 0.30).
#' @param cover_transform "sqrt" (default), "none", or "log1p".
#' @param compare "max" (default), "one_vs_rest", or "pairwise".
#' @return A list with:
#'   - `scores`: matrix "(plots by groups)" of cover-weighted scores,
#'   - `assignments`: data.frame with chosen group per plot and diagnostics,
#'   - `groups`: list of species vectors per group.
#' @export
assign_releves_strict <- function(x, vegmatrix, phi = 0.30,
                                  cover_transform = c("sqrt", "none", "log1p"),
                                  compare = c("max", "one_vs_rest", "pairwise")) {
  cover_transform <- match.arg(cover_transform)
  compare <- match.arg(compare)

  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  spp_all <- colnames(CS)

  # species order identical to plotting & grouping
  ord <- order(H)
  Species.sort <- vapply(seq_len(ncol(CS)), function(i) paste(CS[ord, i], collapse = ""), "")
  species_order <- order(Species.sort)

  # keep only "topmost" clusters whose height >= phi
  idx <- which(H >= phi)
  if (!length(idx)) stop("No clusters at or above phi = ", phi)
  children <- unique(as.integer(CM[idx, ][CM[idx, ] > 0]))
  top_idx <- sort(setdiff(idx, intersect(idx, children)))
  if (!length(top_idx)) stop("No top-level clusters found at phi = ", phi)

  # build strict groups (species sets)
  groups <- lapply(top_idx, function(i) spp_all[CS[i, ] == 1])
  names(groups) <- paste0("G", seq_along(groups))

  # align species between vegmatrix and groups
  vm <- as.matrix(vegmatrix)
  if (is.null(colnames(vm))) stop("vegmatrix must have species column names.")
  common <- intersect(colnames(vm), spp_all)
  if (!length(common)) stop("No overlapping species between vegmatrix and clustering result.")
  vm <- vm[, common, drop = FALSE]

  # weight covers
  W <- .cover_transform(vm, cover_transform)

  # score per plot per group: sum of weights for species in the group
  S <- matrix(0, nrow = nrow(W), ncol = length(groups),
              dimnames = list(rownames(W), names(groups)))
  for (g in seq_along(groups)) {
    sp <- intersect(groups[[g]], colnames(W))
    if (length(sp)) S[, g] <- rowSums(W[, sp, drop = FALSE])
  }

  assignments <- .pick_groups(S, compare = compare)
  list(scores = S, assignments = assignments, groups = groups)
}

# --- 2) Fuzzy assignment using phi weights ---------------------------------

#' Assign releves to fuzzy Cocktail groups (phi-weighted)
#'
#' @param fuzzy_result Output of `cocktail_fuzzy()` (columns: species, group, phi).
#' @param vegmatrix numeric matrix/data.frame: plots x species (covers).
#' @param cover_transform "sqrt" (default), "none", or "log1p".
#' @param weight_by_phi logical; if TRUE (default) multiply cover weights by phi.
#' @param compare "max" (default), "one_vs_rest", or "pairwise".
#' @param phi_min optional filter: only memberships with phi >= phi_min are used
#'   (default NULL = use all).
#' @return A list with:
#'   - `scores`: matrix "(plots x groups)" of phi-cover-weighted scores,
#'   - `assignments`: data.frame with chosen group per plot and diagnostics.
#' @export
assign_releves_fuzzy <- function(fuzzy_result, vegmatrix,
                                 cover_transform = c("sqrt", "none", "log1p"),
                                 weight_by_phi = TRUE,
                                 compare = c("max", "one_vs_rest", "pairwise"),
                                 phi_min = NULL) {
  cover_transform <- match.arg(cover_transform)
  compare <- match.arg(compare)
  stopifnot(all(c("species", "group", "phi") %in% names(fuzzy_result)))

  fr <- fuzzy_result
  if (!is.null(phi_min)) fr <- fr[fr$phi >= phi_min, , drop = FALSE]
  if (!nrow(fr)) stop("No memberships in fuzzy_result after filtering.")

  vm <- as.matrix(vegmatrix)
  if (is.null(colnames(vm))) stop("vegmatrix must have species column names.")

  # keep only species present in data
  fr <- fr[fr$species %in% colnames(vm), , drop = FALSE]
  if (!nrow(fr)) stop("No overlapping species between fuzzy_result and vegmatrix.")

  groups <- sort(unique(fr$group))
  W <- .cover_transform(vm, cover_transform)

  # Build scores: for each group, sum over species (cover_weight * phi or 1)
  S <- matrix(0, nrow = nrow(W), ncol = length(groups),
              dimnames = list(rownames(W), groups))

  # split fr by group and accumulate
  by_g <- split(fr, fr$group)
  for (g in names(by_g)) {
    gg <- by_g[[g]]
    sp <- intersect(gg$species, colnames(W))
    if (!length(sp)) next
    phi_vec <- gg$phi[match(sp, gg$species)]
    if (weight_by_phi) {
      # multiply each species' cover weight by its phi
      S[, g] <- rowSums( sweep(W[, sp, drop = FALSE], 2, phi_vec, `*`) )
    } else {
      # treat membership as unweighted (presence/cover only)
      S[, g] <- rowSums(W[, sp, drop = FALSE])
    }
  }

  assignments <- .pick_groups(S, compare = compare)
  list(scores = S, assignments = assignments)
}
