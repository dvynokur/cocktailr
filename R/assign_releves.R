#' Assign plots to Cocktail groups (strict or fuzzy)
#'
#' @description
#' Assign each plot to one parent Cocktail group at a given \eqn{\phi} cut.
#' In **strict** mode, groups are scored by the summed (optionally transformed)
#' covers of their member species. In **fuzzy** mode, species are weighted by
#' their \eqn{\phi} values against the parent groups (computed from
#' \code{x$Plot.cluster}), and plot scores are cover–\eqn{\phi} weighted sums.
#'
#' @param x A Cocktail result (list with \code{Cluster.species}, \code{Cluster.merged},
#'   \code{Cluster.height}; typically from \code{cocktail_cluster_new()}).
#' @param vegmatrix Numeric matrix/data.frame of covers with plots in rows and species
#'   in columns. Missing values are treated as 0.
#' @param mode Assignment mode: \code{"strict"} or \code{"fuzzy"}.
#' @param phi Numeric \eqn{\phi} cut (between 0 and 1) used to define parent clusters.
#' @param cover_transform Transformation applied to covers before scoring:
#'   \code{"sqrt"} (default), \code{"none"}, or \code{"log1p"}.
#' @param compare Rule for choosing a winner per plot:
#'   \itemize{
#'     \item \code{"max"} — assign only if there is a unique maximum score; otherwise \code{"not assigned"}.
#'     \item \code{"one_vs_rest"} — assign only if the best score is strictly greater than
#'           the sum of all other groups’ scores; otherwise \code{"not assigned"}.
#'   }
#' @param min_cover Numeric threshold for **raw** (pre-transform) cover, expressed in percent
#'   of the plot (between 0 and 100). A plot can be assigned to a group only if the sum of raw
#'   covers of that group’s species in the plot is at least this value. Default is 0.
#' @param min_group_size Integer \eqn{\ge} 1. After initial assignment, any group with fewer
#'   than this many plots is relabeled \code{"not assigned"}. Default is 1 (no collapsing).
#' @param phi_mode (Fuzzy mode only) how to treat per-species \eqn{\phi} weights:
#'   \code{"thresh"} (set \eqn{\phi <} \code{phi} to 0),
#'   \code{"positive"} (set \eqn{\phi \le 0} to 0; default),
#'   or \code{"all"} (keep raw values; negatives penalize).
#'
#' @details
#' Parent clusters are the “topmost” nodes with height \eqn{\ge} \code{phi} (i.e., nodes above the
#' cut that have no ancestor also above the cut).
#' In strict mode, the score for a plot–group pair is the sum of transformed covers across species
#' that belong to that parent group.
#' In fuzzy mode, a species-by-nodes \eqn{\phi} matrix is computed for all internal nodes from
#' \code{x$Plot.cluster}, restricted to the parent groups at \code{phi}, optionally filtered by
#' \code{phi_mode}, and plot scores are obtained by multiplying transformed covers by these \eqn{\phi}
#' weights. The \code{compare} rule then determines whether a unique winner exists for each plot.
#'
#' The \code{min_cover} check is applied in both modes (using raw, untransformed covers); ineligible
#' plot–group pairs are discarded prior to comparison. Finally, \code{min_group_size} can collapse
#' very small winning groups to \code{"not assigned"}.
#'
#' @return A named character vector giving the assigned group for each plot. An attribute
#'   \code{"details"} is attached with inputs and diagnostics, including the score matrix and
#'   a data frame of per-plot summaries.
#'
#' @seealso \code{\link{cocktail_fuzzy}}
#' @export

assign_releves <- function(
    x,
    vegmatrix,
    mode            = c("strict", "fuzzy"),
    phi             = 0.30,
    cover_transform = c("sqrt", "none", "log1p"),
    compare         = c("max", "one_vs_rest"),
    min_cover       = 0,
    min_group_size  = 1,
    phi_mode        = NULL   # fuzzy-only; if provided in strict → error
) {
  ## ---- helpers ----
  .cover_transform <- function(xx, how) {
    if (how == "sqrt") sqrt(pmax(xx, 0))
    else if (how == "log1p") log1p(pmax(xx, 0))
    else xx
  }
  .top_parent_clusters <- function(CM, H, cut) {
    idx <- which(H >= cut); if (!length(idx)) return(integer(0))
    kids <- CM[idx, , drop = FALSE]
    children <- unique(as.integer(kids[kids > 0]))
    sort(setdiff(idx, intersect(idx, children)))
  }
  .is_cocktail_object <- function(obj) {
    is.list(obj) && all(c("Cluster.species","Cluster.merged","Cluster.height") %in% names(obj))
  }

  ## ---- args / checks ----
  mode            <- match.arg(mode)
  cover_transform <- match.arg(cover_transform)
  compare         <- match.arg(compare)

  if (!.is_cocktail_object(x)) {
    stop("x must be a Cocktail object (with Cluster.species, Cluster.merged, Cluster.height).")
  }
  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  spp_all <- colnames(CS)

  # vegmatrix
  vm <- as.matrix(vegmatrix)
  vm[is.na(vm)] <- 0
  storage.mode(vm) <- "double"
  if (is.null(colnames(vm))) stop("vegmatrix must have species column names.")
  if (is.null(rownames(vm))) rownames(vm) <- paste0("plot_", seq_len(nrow(vm)))

  # clip thresholds
  min_cover <- max(0, min(100, as.numeric(min_cover)))
  min_group_size <- as.integer(min_group_size)
  if (is.na(min_group_size) || min_group_size < 1) min_group_size <- 1L

  # strict must not get phi_mode
  if (mode == "strict" && !is.null(phi_mode)) {
    stop("`phi_mode` is only allowed when mode = 'fuzzy'.")
  }

  # parent clusters at/above phi
  top_idx <- .top_parent_clusters(CM, H, phi)
  if (!length(top_idx)) stop("No parent clusters at or above phi = ", phi)
  group_names <- paste0("c_", top_idx)

  # align species
  common <- intersect(colnames(vm), spp_all)
  if (!length(common)) stop("No overlapping species between vegmatrix and clustering result.")
  vm <- vm[, common, drop = FALSE]
  W  <- .cover_transform(vm, cover_transform)  # transformed covers
  vm_raw <- vm                                 # raw covers for min_cover

  # group species sets (used for scoring in strict and for min_cover in both)
  groups_species <- lapply(top_idx, function(i) spp_all[CS[i, ] == 1L])
  names(groups_species) <- group_names

  # precompute RAW cover per group (min_cover eligibility)
  S_raw <- matrix(0, nrow = nrow(vm_raw), ncol = length(groups_species),
                  dimnames = list(rownames(vm_raw), group_names))
  for (g in seq_along(groups_species)) {
    sp <- intersect(groups_species[[g]], colnames(vm_raw))
    if (length(sp)) S_raw[, g] <- rowSums(vm_raw[, sp, drop = FALSE])
  }
  S_raw <- pmin(S_raw, 100)  # cap to 100%

  ## ---- scores ----
  if (mode == "strict") {
    # sum of transformed covers over group species
    S <- matrix(0, nrow = nrow(W), ncol = length(groups_species),
                dimnames = list(rownames(W), group_names))
    for (g in seq_along(groups_species)) {
      sp <- intersect(groups_species[[g]], colnames(W))
      if (length(sp)) S[, g] <- rowSums(W[, sp, drop = FALSE])
    }
  } else {
    # fuzzy: compute phi from Cocktail core (Plot.cluster) for ALL nodes, then subset to parents
    if (is.null(x$Plot.cluster)) stop("x$Plot.cluster is required for fuzzy mode.")

    # presence/absence (0/1 doubles)
    X <- (vm > 0); storage.mode(X) <- "double"     # N × nsp
    N <- nrow(X); nsp <- ncol(X)

    G <- as.matrix(x$Plot.cluster > 0); storage.mode(G) <- "double"  # N × n_nodes
    if (nrow(G) != N) stop("Internal mismatch: Plot.cluster has different number of plots.")
    n_nodes <- ncol(G)

    # a = t(X) %*% G  (nsp × n_nodes)
    a  <- crossprod(X, G)
    p  <- colSums(X)     # nsp
    g1 <- colSums(G)     # n_nodes

    # broadcast margins
    b <- matrix(p,  nrow = nsp, ncol = n_nodes) - a
    c <- matrix(g1, nrow = nsp, ncol = n_nodes, byrow = TRUE) - a
    d <- (N - matrix(g1, nrow = nsp, ncol = n_nodes, byrow = TRUE)) - b

    den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
    Phi_all <- (a * d - b * c) / den
    Phi_all[!is.finite(den) | den <= 0] <- 0

    colnames(Phi_all) <- paste0("c_", seq_len(n_nodes))
    rownames(Phi_all) <- colnames(X)

    # keep only parent columns at the cut
    keep_cols <- intersect(colnames(Phi_all), group_names)
    Phi <- Phi_all[, keep_cols, drop = FALSE]

    # phi_mode default + filtering
    phi_mode <- if (is.null(phi_mode)) "positive" else match.arg(phi_mode, c("thresh","positive","all"))
    if (phi_mode == "thresh") {
      Phi[Phi <  phi] <- 0
    } else if (phi_mode == "positive") {
      Phi[Phi <= 0] <- 0
    }
    # scores: cover weights × phi
    S <- as.matrix(W %*% Phi)
    # ensure final column order matches group_names
    S <- S[, group_names, drop = FALSE]
  }

  ## ---- choose winner per plot ----
  # apply min_cover eligibility: ineligible groups get -Inf
  if (min_cover > 0) {
    inelig <- S_raw < min_cover
    S[inelig] <- -Inf
  }

  rn <- rownames(S); cn <- colnames(S)
  best_idx   <- max.col(S, ties.method = "first")
  best_score <- S[cbind(seq_len(nrow(S)), best_idx)]
  best_name  <- cn[best_idx]

  # second best (diagnostics)
  second_score <- apply(S, 1, function(v) {
    if (length(v) == 1) 0 else sort(v, decreasing = TRUE)[2L]
  })

  assigned <- rep("not assigned", length(best_idx))

  if (compare == "max") {
    ties_per_row <- apply(S, 1, function(v) sum(v == max(v)))
    is_unique <- (ties_per_row == 1) & is.finite(best_score)
    assigned[is_unique] <- best_name[is_unique]
  } else { # "one_vs_rest"
    row_sums <- rowSums(S)
    rest     <- pmax(row_sums - best_score, 0)
    ok       <- is.finite(best_score) & (best_score > rest)
    assigned[ok] <- best_name[ok]
  }

  # collapse tiny groups
  collapsed <- character(0)
  if (min_group_size > 1) {
    counts <- table(assigned)
    small  <- names(counts)[counts < min_group_size & names(counts) != "not assigned"]
    if (length(small)) {
      assigned[assigned %in% small] <- "not assigned"
      collapsed <- small
    }
  }

  names(assigned) <- rn
  assignments <- data.frame(
    plot           = rn,
    assigned_group = assigned,
    best_score     = best_score,
    second_best    = second_score,
    margin         = best_score - second_score,
    stringsAsFactors = FALSE
  )

  out <- stats::setNames(assigned, rn)
  attr(out, "details") <- list(
    mode             = mode,
    phi              = phi,
    compare          = compare,
    cover_transform  = cover_transform,
    min_cover        = min_cover,
    min_group_size   = min_group_size,
    collapsed_groups = collapsed,
    scores           = S,
    assignments      = assignments
  )
  out
}
