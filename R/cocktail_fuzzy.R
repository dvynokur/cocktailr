#' Fuzzy Cocktail memberships (species vs all internal nodes)
#'
#' @description
#' Compute a species-by-nodes matrix of \eqn{\phi} association coefficients for
#' **all** internal nodes (clusters) in a fitted Cocktail tree. Group membership
#' per node is taken from `x$Plot.cluster`, while species occurrences come
#' from `vegmatrix` after binarization (> 0 becomes 1; `NA` is 0).
#'
#' @param x Cocktail result (must contain `Cluster.species`, `Cluster.merged`,
#'   `Cluster.height`, `Plot.cluster`, and carry `species` and `plots` names).
#' @param vegmatrix A plots-by-species numeric matrix used at clustering time
#'   (or the same grid). Values > 0 are treated as presence; `NA` is set to 0.
#'   Row names (plots) and column names (species) are used to align with `x`.
#'
#' @details
#' - The function aligns rows/columns to the plots and species stored in `x`.
#'   Extra plots in `vegmatrix` are ignored; extra species with any nonzero
#'   values cause an error (to avoid mismatched inputs).
#' - The \eqn{\phi} for every (species, node) pair is computed from the
#'   2×2 table formed by species presence and node membership across plots.
#'   Invalid or zero denominators yield \eqn{\phi}=0.
#'
#' @return A numeric matrix of size (species × nodes) with columns named
#'   `"c_<node_id>"`.
#'
#' @seealso [assign_releves()], [species_in_clusters()]
#' @export

cocktail_fuzzy <- function(x, vegmatrix) {
  # ---- sanity on x ----
  need <- c("Cluster.species", "Cluster.merged", "Cluster.height", "Plot.cluster")
  miss <- setdiff(need, names(x))
  if (length(miss)) stop("x is missing: ", paste(miss, collapse = ", "))

  CS <- x$Cluster.species
  PC <- x$Plot.cluster
  sp_x <- if (!is.null(x$species)) x$species else colnames(CS)
  pl_x <- if (!is.null(x$plots))   x$plots   else rownames(PC)
  if (is.null(sp_x) || is.null(pl_x))
    stop("x must carry `species` and `plots` names (set by cocktail_cluster_*()).")

  # ---- ingest & binarize vegmatrix (double-safe) ----
  vm <- as.matrix(vegmatrix)
  if (is.null(colnames(vm)) || is.null(rownames(vm)))
    stop("vegmatrix must have species (colnames) and plot (rownames).")
  vm[is.na(vm)] <- 0
  storage.mode(vm) <- "double"

  # align to clustering grid
  if (!all(pl_x %in% rownames(vm))) {
    missing_pl <- setdiff(pl_x, rownames(vm))
    stop("vegmatrix is missing plots used in clustering: ",
         paste(head(missing_pl, 10), collapse = ", "),
         if (length(missing_pl) > 10) " ..." else "")
  }
  if (!all(sp_x %in% colnames(vm))) {
    missing_sp <- setdiff(sp_x, colnames(vm))
    stop("vegmatrix is missing species used in clustering: ",
         paste(head(missing_sp, 10), collapse = ", "),
         if (length(missing_sp) > 10) " ..." else "")
  }
  vm <- vm[pl_x, sp_x, drop = FALSE]

  # presence/absence (N × nsp), as numeric 0/1
  X <- (vm > 0)
  storage.mode(X) <- "double"

  N      <- nrow(X)
  nsp    <- ncol(X)
  n_nodes <- nrow(CS)             # internal nodes = n - 1

  # ---- group membership for ALL nodes from Cocktail core ----
  G  <- as.matrix(PC > 0)         # N × n_nodes (logical)
  storage.mode(G) <- "double"
  if (nrow(G) != N || ncol(G) != n_nodes)
    stop("Internal mismatch: Plot.cluster size differs from data.")

  # ---- co-occurrence counts for ALL (species,node) at once ----
  # a = t(X) %*% G  → nsp × n_nodes
  a <- crossprod(X, G)

  # margins
  p  <- colSums(X)     # length nsp (species totals = a+b)
  g1 <- colSums(G)     # length n_nodes (group sizes)

  # broadcasted matrices
  b <- matrix(p,  nrow = nsp, ncol = n_nodes) - a
  c <- matrix(g1, nrow = nsp, ncol = n_nodes, byrow = TRUE) - a
  d <- (N - matrix(g1, nrow = nsp, ncol = n_nodes, byrow = TRUE)) - b

  den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
  phi <- (a * d - b * c) / den
  phi[!is.finite(den) | den <= 0] <- 0

  colnames(phi) <- paste0("c_", seq_len(n_nodes))
  rownames(phi) <- sp_x

  GI <- data.frame(
    col        = colnames(phi),
    cluster_id = seq_len(n_nodes),
    n_plots    = as.numeric(g1),
    stringsAsFactors = FALSE
  )
  attr(phi, "group_info") <- GI

  phi
}
