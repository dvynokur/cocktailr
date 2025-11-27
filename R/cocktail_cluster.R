#' Cocktail clustering (sparse matrix)
#'
#' Fast Cocktail agglomeration for a **plots × species** table.
#'
#' @description
#' This implementation:
#' - **Binarizes** the input for clustering: values > 0 become 1; values ≤ 0 or `NA` become 0.
#' - **Drops empty species** (all-zero columns) before clustering (and keeps `vegmatrix` aligned).
#' - Computes the association coefficient (“phi”) each round from one sparse
#'   crossproduct for speed and exactness.
#' - Uses a **fixed, reproducible tie order**: when several pairs share the same
#'   maximum phi at a step, they are processed in the same order that R fills the
#'   lower-triangular distance matrix (scan by increasing column, then row).
#' - Stores `Plot.cluster` as a **sparse `dgCMatrix`** (from the **Matrix** package),
#'   containing either binary membership or relative cover per plot and node.
#'   If you need a base R matrix, convert manually, e.g.:
#'   `plot_cluster_dense <- as.matrix(res$Plot.cluster)`.
#' - Optionally writes **relative cover** into `Plot.cluster` instead of binary membership
#'   via `plot_values = "rel_cover"`; relative cover is defined as
#'   (sum of cluster species covers per plot) / (total cover of the plot),
#'   and values are zeroed for plots not meeting the m-threshold (cluster membership).
#' - Optionally computes **species–cluster association coefficients** (`Species.cluster.phi`),
#'   a species × nodes matrix of \eqn{\phi} between species presence and node membership,
#'   using sparse crossproducts internally.
#'
#' @param vegmatrix A matrix or data frame with **plots in rows** and **species in columns**.
#'   This is used twice:
#'   (1) it is **binarized** internally to drive the clustering, and
#'   (2) its **original numeric values** (with `NA` treated as 0) are used to compute
#'       relative cover when `plot_values = "rel_cover"`.
#' @param progress  Logical; show a text progress bar (default `TRUE`).
#' @param plot_values Character; one of `c("binary", "rel_cover")`.
#'   - `"binary"` (default): `Plot.cluster` stores 0/1 plot membership per merge.
#'   - `"rel_cover"`: `Plot.cluster` stores the **relative cover** per plot and merge:
#'       sum of covers over the cluster’s species divided by the total cover of the plot,
#'       but **only** for plots that meet the current merge’s m-threshold (membership);
#'       non-member plots or plots with zero total cover are set to 0.
#' @param species_cluster_phi Logical; if `TRUE`, compute and return
#'   `Species.cluster.phi`, a species × nodes matrix of \eqn{\phi} association
#'   coefficients between species presence (from binarized `vegmatrix`) and node
#'   membership (from `Plot.cluster > 0`). Default `FALSE` to avoid the extra
#'   computation and memory cost.
#'
#' @return
#' A list of class `"cocktail"` with:
#' \itemize{
#'   \item `Cluster.species`       — integer matrix (n_merges × n_species): species membership per merge.
#'   \item `Cluster.info`          — integer matrix (n_merges × 2): columns `k` (cluster size) and `m` (threshold).
#'   \item `Plot.cluster`          — **`dgCMatrix`** (n_plots × n_merges): plot values per merge
#'                                   (0/1 for `"binary"`, relative cover for `"rel_cover"`).
#'   \item `Cluster.merged`        — integer matrix (n_merges × 2): left/right children per merge
#'                                   (negative = original species index; positive = earlier merge index).
#'   \item `Cluster.height`        — numeric vector length n_merges: phi at each merge.
#'   \item `Species.cluster.phi`   — (optional) numeric matrix (species × nodes) of \eqn{\phi} association
#'                                   coefficients between each species and each internal node
#'                                   (columns named `"c_<node_id>"`), with an attribute
#'                                   `"group_info"` giving node sizes. `NULL` if
#'                                   `species_cluster_phi = FALSE`.
#'   \item `species`               — character vector of species names kept after cleaning.
#'   \item `plots`                 — character vector of plot names.
#' }
#'
#' @details
#' - Binarization and removal of empty species happen internally and only affect the
#'   set of columns that contribute to clustering. All returned components are aligned
#'   to the species that had at least one presence after cleaning.
#' - For `plot_values = "rel_cover"`, relative cover is computed from the **original**
#'   `vegmatrix` values (after converting `NA` to 0). For each merge, the function:
#'   (1) identifies the cluster’s species, (2) sums their covers per plot,
#'   (3) divides by the total cover in that plot (sum over all species), and
#'   (4) **zeroes** values for plots that do not meet the m-threshold
#'       (i.e., are not assigned to that merge) or have zero total cover.
#' - Basic checks on the cover scale:
#'   if input values appear **binary** (only 0/1) or contain **non-numeric codes**
#'   (e.g., `+, r, 2a`), the function warns and **falls back to `"binary"`** output.
#'   If values look **ordinal** (e.g., 1..6 / 1..10), the function warns but proceeds
#'   to compute relative cover, noting that percentage cover is recommended.
#' - `Species.cluster.phi` is computed from a 2×2 table for every (species,node) pair,
#'   using species presence (from binarized `vegmatrix`) and node membership
#'   (from `Plot.cluster > 0`). Sparse crossproducts (`Matrix::crossprod`) are used
#'   to obtain co-occurrence counts efficiently; invalid or zero denominators yield
#'   \eqn{\phi}=0.
#'
#' @import Matrix
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
cocktail_cluster <- function(
    vegmatrix,
    progress = TRUE,
    plot_values = c("binary", "rel_cover"),
    species_cluster_phi = FALSE
) {
  plot_values <- match.arg(plot_values)

  ## ---- input checks & setup ----
  if (!is.matrix(vegmatrix) && !is.data.frame(vegmatrix)) {
    stop("vegmatrix must be a matrix or data.frame with plots in rows and species in columns.")
  }

  # Keep original covers for relative cover output and species-cluster phi
  vm_raw <- as.matrix(vegmatrix)
  vm_raw[is.na(vm_raw)] <- 0

  # Detect cover scale for warnings/forcing behavior
  detect_cover_scale <- function(M) {
    vals <- unique(as.vector(M)); vals <- vals[!is.na(vals)]
    if (!length(vals)) return(list(type="unknown", note=NULL))
    num_try <- suppressWarnings(as.numeric(vals))
    non_num <- is.na(num_try)
    if (any(non_num)) return(list(type="non_numeric", note="Non-numeric cover codes detected (e.g., '+', 'r', '2a')."))
    u <- sort(unique(num_try))
    if (all(u %in% c(0,1))) return(list(type="binary", note="Cover data appear to be binary (0/1)."))
    all_int <- all(abs(u - round(u)) < .Machine$double.eps^0.5)
    if (all_int && length(u) <= 10) return(list(type="ordinal", note="Cover data appear ordinal (small integer scale)."))
    list(type="numeric", note=NULL)
  }
  scale_info <- detect_cover_scale(vm_raw)

  # If user asked for relative cover but data unsuitable, warn and force binary
  if (plot_values != "binary") {
    if (scale_info$type %in% c("binary","non_numeric")) {
      warning(sprintf(
        "%s Using binary Plot.cluster instead.",
        if (scale_info$type == "binary") {
          "Cover data are binary (0/1); relative covers are not meaningful."
        } else {
          "Non-numeric cover codes detected; relative covers not computed."
        }
      ))
      plot_values <- "binary"
    } else if (scale_info$type == "ordinal") {
      warning("Cover data look ordinal (e.g., 1..6 / 1..10). Proceeding with relative cover, but percentage cover is recommended.")
    }
  }

  # Binarize for clustering
  vm <- as.matrix(vegmatrix)
  vm[is.na(vm)] <- 0
  vm <- (vm > 0) * 1L

  plots   <- rownames(vm); if (is.null(plots))   plots   <- as.character(seq_len(nrow(vm)))
  species <- colnames(vm); if (is.null(species)) species <- as.character(seq_len(ncol(vm)))

  # Drop empty species columns (and align vm_raw)
  keep <- colSums(vm) > 0L
  if (!all(keep)) {
    vm      <- vm[, keep, drop = FALSE]
    vm_raw  <- vm_raw[, keep, drop = FALSE]
    species <- species[keep]
  }

  N <- nrow(vm); n <- ncol(vm)
  if (n < 2L) stop("After dropping empty species, need at least 2 species (columns).")
  if (N < 1L) stop("Need at least 1 plot (row).")

  # species-only sparse matrix
  X0 <- Matrix::Matrix(vm, sparse = TRUE)  # dgCMatrix (0/1)
  rm(vm)

  # global frequencies for Expected.plot.freq (from original species)
  p.freq <- as.numeric(Matrix::colSums(X0)) / N
  q.freq <- 1 - p.freq

  # Total cover per plot (denominator for relative cover)
  plot_totals <- rowSums(vm_raw, na.rm = TRUE)

  ## ---- helpers -------------------------------------------------------------
  Expected.plot.freq_ <- function(species.in.cluster) {
    K <- length(species.in.cluster)
    Exp <- array(0, K + 1L)
    Exp_inter <- array(0, K + 1L)
    Exp_inter[1L] <- 1
    for (j in 1L:K) {
      s <- species.in.cluster[j]
      Exp[1L] <- Exp_inter[1L] * q.freq[s]
      if (j > 1L) {
        for (k in 1L:(j - 1L)) {
          Exp[k + 1L] <- Exp_inter[k]     * p.freq[s] +
            Exp_inter[k + 1L] * q.freq[s]
        }
      }
      Exp[j + 1L] <- Exp_inter[j] * p.freq[s]
      for (k in 1L:(j + 1L)) Exp_inter[k] <- Exp[k]
    }
    Exp
  }

  Compare.obs.exp.freq_ <- function(Obs.freq, Exp.freq) {
    Obs <- if (is.matrix(Obs.freq)) as.vector(Obs.freq[, 1L]) else as.vector(Obs.freq)
    K   <- length(Obs) - 1L
    Cum.obs <- array(0, K + 1L)
    Cum.exp <- array(0, K + 1L)
    Cum.obs[K + 1L] <- Obs[K + 1L]
    Cum.exp[K + 1L] <- Exp.freq[K + 1L]
    m <- 1L
    m.found <- -1L
    if (Cum.obs[K + 1L] > Cum.exp[K + 1L]) { m <- K; m.found <- 0L }
    for (j in K:1L) {
      Cum.obs[j] <- Cum.obs[j + 1L] + Obs[j]
      if (m.found == -1L && Cum.obs[j] > 0) m.found <- 0L
      Cum.exp[j] <- Cum.exp[j + 1L] + Exp.freq[j]
      if (j > 1L && m.found == 0L && Cum.exp[j] > Cum.obs[j]) { m <- j; m.found <- 1L }
    }
    m
  }

  .lower_tri_index_vec <- function(i, j, m) {
    x <- i - 1L
    (x * m) - (x * (x + 1L)) / 2 + (j - i)
  }

  # phi via one sparse crossproduct
  phi_max_pairs_crossprod_distorder_ <- function(X) {
    stopifnot(inherits(X, "dgCMatrix"))
    m <- ncol(X)
    if (m < 2L) return(list(max_phi = 0, pairs = matrix(integer(0), ncol = 2)))

    Nn <- nrow(X)
    p  <- as.numeric(Matrix::colSums(X))  # counts

    A <- Matrix::crossprod(X)             # dsCMatrix symmetric, upper triangle stored
    Matrix::diag(A) <- 0L
    if (length(A@x) == 0L) {
      # no co-occurrences → need fallback to include a==0 pairs
      return(list(max_phi = -Inf, pairs = matrix(integer(0), ncol = 2)))
    }

    # indices aligned with A@x (upper triangle, compressed-by-column storage)
    j_idx <- rep.int(seq_len(m), diff(A@p))  # column index per nonzero
    i_idx <- A@i + 1L                        # row index per nonzero (i < j)
    a     <- A@x                             # co-occurrence counts

    b <- p[i_idx] - a
    c <- p[j_idx] - a
    d <- Nn - a - b - c

    den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
    phi <- ifelse(den > 0, (a * d - b * c) / den, 0)
    phi[!is.finite(phi)] <- 0

    max_phi <- max(phi)
    keep <- which(phi == max_phi)
    if (!length(keep)) return(list(max_phi = -Inf, pairs = matrix(integer(0), ncol = 2)))

    # sort ties by the original lower-triangle linear index (exact which()-order)
    idx <- .lower_tri_index_vec(i_idx[keep], j_idx[keep], m)
    ord <- order(idx)
    pairs <- cbind(e1 = i_idx[keep][ord], e2 = j_idx[keep][ord])
    list(max_phi = max_phi, pairs = pairs)
  }

  # Fallback φ over all pairs, ties sorted by the same lower-triangle index
  phi_max_pairs_fallback_distorder_ <- function(X) {
    stopifnot(inherits(X, "dgCMatrix"))
    m <- ncol(X)
    if (m < 2L) return(list(max_phi = 0, pairs = matrix(integer(0), ncol = 2)))

    Nn <- nrow(X)
    p  <- as.numeric(Matrix::colSums(X))

    best <- -Inf
    out_i <- integer(0)
    out_j <- integer(0)

    for (j in 2L:m) {
      a <- as.numeric(Matrix::t(X[, 1L:(j - 1L), drop = FALSE]) %*% X[, j, drop = FALSE])
      b <- p[1L:(j - 1L)] - a
      c <- p[j]           - a
      d <- Nn - a - b - c
      den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
      phi <- ifelse(den > 0, (a * d - b * c) / den, 0)
      phi[!is.finite(phi)] <- 0

      cur <- max(phi)
      if (cur > best) {
        keep <- which(phi == cur)
        best <- cur
        out_i <- keep
        out_j <- rep.int(j, length(keep))
      } else if (cur == best) {
        keep <- which(phi == best)
        if (length(keep)) {
          out_i <- c(out_i, keep)
          out_j <- c(out_j, rep.int(j, length(keep)))
        }
      }
    }

    if (!length(out_i)) {
      return(list(max_phi = ifelse(is.finite(best), best, 0),
                  pairs   = matrix(integer(0), ncol = 2)))
    }
    idx <- .lower_tri_index_vec(out_i, out_j, m)
    ord <- order(idx)
    pairs <- cbind(e1 = out_i[ord], e2 = out_j[ord])
    list(max_phi = ifelse(is.finite(best), best, 0), pairs = pairs)
  }

  ## ---- outputs ----
  Cluster.species <- matrix(0L, n - 1L, n, dimnames = list(NULL, species))
  Cluster.info    <- matrix(0L, n - 1L, 2L,
                            dimnames = list(as.character(seq_len(n - 1L)), c("k","m")))
  # Plot.cluster as sparse dgCMatrix
  Plot.cluster    <- Matrix::Matrix(0, nrow = N, ncol = n - 1L,
                                    sparse = TRUE,
                                    dimnames = list(plots, NULL))
  Cluster.merged  <- matrix(0L, n - 1L, 2L)
  Cluster.height  <- array(0, n - 1L)

  ## ---- working state ----
  X         <- X0
  col_names <- colnames(X0); if (is.null(col_names)) col_names <- species
  i <- 0L
  name_last_cluster <- NULL

  pb <- if (isTRUE(progress)) utils::txtProgressBar(min = 0, max = n - 1L, style = 3) else NULL
  on.exit({ if (!is.null(pb)) close(pb) }, add = TRUE)

  ## ---- agglomeration loop ----
  while (i <= (n - 2L)) {

    # fast φ; if needed, fallback (ensures a==0 pairs included when max ≤ 0)
    maxres <- phi_max_pairs_crossprod_distorder_(X)
    if (nrow(maxres$pairs) == 0L || !(is.finite(maxres$max_phi) && maxres$max_phi > 0)) {
      maxres <- phi_max_pairs_fallback_distorder_(X)
      if (nrow(maxres$pairs) == 0L) break
    }

    e1 <- maxres$pairs[, 1L]
    e2 <- maxres$pairs[, 2L]
    multiple.max <- length(e1)

    # circularity filter
    compare1 <- as.vector(t(cbind(e1, e2)))
    if (anyDuplicated(compare1) > 2L) {
      keep <- rep(TRUE, multiple.max)
      for (jj in 2L:multiple.max) {
        compare2 <- compare1[1L:(2L * (jj - 1L))]
        k1 <- sum(!is.na(match(compare2, e1[jj])))
        k2 <- sum(!is.na(match(compare2, e2[jj])))
        if (k1 > 0L & k2 > 0L) keep[jj] <- FALSE
      }
      e1 <- e1[keep]; e2 <- e2[keep]
      multiple.max <- length(e1)
      if (multiple.max == 0L) next
    }

    # respect the (n-1) bound
    remaining_merges <- (n - 1L) - i
    if (remaining_merges <= 0L) break
    if (multiple.max > remaining_merges) {
      e1 <- e1[seq_len(remaining_merges)]
      e2 <- e2[seq_len(remaining_merges)]
      multiple.max <- remaining_merges
    }

    i1 <- i + 1L

    # ----- rename only within current endpoints, and evolve within round -----
    end_idx   <- c(e1, e2)
    end_names <- col_names[end_idx]

    for (jj in seq_len(multiple.max)) {
      if (i >= (n - 1L)) break
      i <- i + 1L
      Cluster.height[i] <- maxres$max_phi

      pos_left  <- match(e1[jj], end_idx)
      pos_right <- match(e2[jj], end_idx)

      left_name  <- end_names[pos_left]
      right_name <- end_names[pos_right]

      if (startsWith(left_name, "c_")) {
        cl1 <- as.integer(sub("c_", "", left_name))
        Cluster.merged[i, 1L] <- cl1
        Cluster.species[i, Cluster.species[cl1, ] == 1L] <- 1L
      } else {
        f1 <- match(left_name, species)
        Cluster.merged[i, 1L] <- -f1
        Cluster.species[i, f1] <- 1L
      }

      if (startsWith(right_name, "c_")) {
        cl2 <- as.integer(sub("c_", "", right_name))
        Cluster.merged[i, 2L] <- cl2
        Cluster.species[i, Cluster.species[cl2, ] == 1L] <- 1L
      } else {
        f2 <- match(right_name, species)
        Cluster.merged[i, 2L] <- -f2
        Cluster.species[i, f2] <- 1L
      }

      # evolve endpoint names inside the current batch
      newc <- paste0("c_", i)
      tmp <- end_names
      tmp[tmp == left_name]  <- newc
      tmp[tmp == right_name] <- newc
      end_names <- tmp

      # k, m, plot assignment (computed from original X0)
      k <- sum(Cluster.species[i, ])
      Cluster.info[i, 1L] <- k

      spp_idx <- which(Cluster.species[i, ] > 0L)
      species_in_plot <- as.integer(Matrix::rowSums(X0[, spp_idx, drop = FALSE]))
      Obs_plot_freq <- tabulate(species_in_plot + 1L, nbins = k + 1L)
      Exp_plot_freq <- Expected.plot.freq_(spp_idx) * N
      Cluster.info[i, 2L] <- Compare.obs.exp.freq_(Obs_plot_freq, Exp_plot_freq)

      # Fill Plot.cluster per option
      if (plot_values == "binary") {
        # 0/1 membership
        Plot.cluster[species_in_plot >= Cluster.info[i, 2L], i] <- 1

      } else if (plot_values == "rel_cover") {
        # Relative cover: sum of cluster covers / total cover per plot,
        # zeroed outside membership or when total cover is zero
        cov_block <- vm_raw[, spp_idx, drop = FALSE]
        sums <- rowSums(cov_block, na.rm = TRUE)
        rel <- ifelse(plot_totals > 0, sums / plot_totals, 0)
        rel[species_in_plot < Cluster.info[i, 2L]] <- 0
        Plot.cluster[, i] <- rel
      }
    }

    col_names[end_idx] <- end_names

    n2 <- ncol(X)
    if (n2 == 3L) {
      name_last_cluster <- col_names[-c(e1[1], e2[1])]
    }

    idx_e1_unique <- !duplicated(col_names[e1], fromLast = TRUE)
    idx_e2_unique <- !duplicated(col_names[e2], fromLast = TRUE)
    index.e <- seq.int(i1, i)[idx_e1_unique & idx_e2_unique]

    drop_cols <- sort(unique(c(e1, e2)))
    if (length(index.e)) {
      newcols <- do.call(cbind, lapply(index.e, function(kid) {
        Matrix::Matrix(Plot.cluster[, kid] > 0, sparse = TRUE)
      }))
      colnames(newcols) <- paste0("c_", index.e)
      X <- cbind(X[, -drop_cols, drop = FALSE], newcols)
      col_names <- c(col_names[-drop_cols], colnames(newcols))
    } else {
      X <- X[, -drop_cols, drop = FALSE]
      col_names <- col_names[-drop_cols]
    }

    if (ncol(X) == 2L && length(name_last_cluster) == 1L) {
      col_names[1L] <- name_last_cluster
    }

    # φ = 0 tail: once every plot is in the current cluster, finish remaining merges at height 0
    if (sum(Plot.cluster[, i] > 0) == N) {
      remaining <- (n - 1L) - i
      if (remaining > 0L) {
        for (j in (i + seq_len(remaining))) {
          # sanity: j within 1..n-1
          stopifnot(j >= 1L, j <= nrow(Cluster.merged))

          g1  <- ncol(X)
          cl1 <- as.integer(sub("c_", "", col_names[g1]))
          Cluster.merged[j, 1L] <- cl1

          g2  <- 1L
          cl2 <- as.integer(sub("c_", "", col_names[g2]))
          Cluster.merged[j, 2L] <- cl2

          Cluster.height[j] <- 0
          Cluster.info[j, 1L] <- sum(Cluster.species[j, ])
          Cluster.info[j, 2L] <- 1L

          Cluster.species[j, Cluster.species[cl1, ] == 1L] <- 1L
          Cluster.species[j, Cluster.species[cl2, ] == 1L] <- 1L

          # Tail: write Plot.cluster
          if (plot_values == "binary") {

            Plot.cluster[, j] <- 1

          } else if (plot_values == "rel_cover") {

            spp_idx <- which(Cluster.species[j, ] > 0L)
            cov_block <- vm_raw[, spp_idx, drop = FALSE]
            sums <- if (length(spp_idx)) rowSums(cov_block, na.rm = TRUE) else rep(0, N)
            species_in_plot_tail <- as.integer(Matrix::rowSums(X0[, spp_idx, drop = FALSE]))
            rel <- ifelse(plot_totals > 0, sums / plot_totals, 0)
            rel[species_in_plot_tail < 1L] <- 0
            Plot.cluster[, j] <- rel
          }

          X <- cbind(X, Matrix::Matrix(Plot.cluster[, j] > 0, sparse = TRUE))
          colnames(X)[ncol(X)] <- paste0("c_", j)
          X <- X[, -c(g1, g2), drop = FALSE]
          col_names <- c(col_names, paste0("c_", j))
          col_names <- col_names[-c(g1, g2)]
        }
        # advance i by exactly how many merges we just filled
        i <- i + remaining
      }
    }

    if (!is.null(pb)) utils::setTxtProgressBar(pb, min(i, n - 1L))
    if (i >= (n - 1L)) break
  }

  ## ---- species × node phi matrix (optional species–cluster associations) ----
  Species.cluster.phi <- NULL
  if (isTRUE(species_cluster_phi)) {
    n_nodes <- nrow(Cluster.species)

    # species presence (plots × species, 0/1) as numeric sparse (dgCMatrix)
    X_f <- Matrix::Matrix((vm_raw > 0) * 1, sparse = TRUE)
    # node membership (plots × nodes, 0/1) as numeric sparse (dgCMatrix)
    G_f <- Matrix::Matrix((Plot.cluster > 0) * 1, sparse = TRUE)

    if (nrow(X_f) != nrow(G_f) || ncol(G_f) != n_nodes) {
      stop("Internal mismatch computing Species.cluster.phi: dimensions do not align.")
    }

    # co-occurrences via sparse crossprod
    a_sc <- Matrix::crossprod(X_f, G_f)   # nsp × n_nodes (Matrix)
    a_sc <- as.matrix(a_sc)               # dense numeric for phi algebra

    # margins (sparse-aware)
    p_sc  <- Matrix::colSums(X_f)         # species totals
    g1_sc <- Matrix::colSums(G_f)         # node sizes

    nsp    <- length(p_sc)
    Nf     <- nrow(X_f)
    n_nodes_check <- length(g1_sc)
    stopifnot(nsp == nrow(a_sc), n_nodes_check == ncol(a_sc))

    # broadcasted
    b_sc <- matrix(p_sc,  nrow = nsp, ncol = n_nodes_check) - a_sc
    c_sc <- matrix(g1_sc, nrow = nsp, ncol = n_nodes_check, byrow = TRUE) - a_sc
    d_sc <- (Nf - matrix(g1_sc, nrow = nsp, ncol = n_nodes_check, byrow = TRUE)) - b_sc

    den_sc <- sqrt((a_sc + c_sc) * (b_sc + d_sc) * (a_sc + b_sc) * (c_sc + d_sc))
    phi_sc <- (a_sc * d_sc - b_sc * c_sc) / den_sc
    phi_sc[!is.finite(den_sc) | den_sc <= 0] <- 0

    colnames(phi_sc) <- paste0("c_", seq_len(n_nodes_check))
    rownames(phi_sc) <- species

    GI <- data.frame(
      col        = colnames(phi_sc),
      cluster_id = seq_len(n_nodes_check),
      n_plots    = as.numeric(g1_sc),
      stringsAsFactors = FALSE
    )
    attr(phi_sc, "group_info") <- GI

    Species.cluster.phi <- phi_sc
  }

  ## ---- result ----
  res <- list(
    Cluster.species      = Cluster.species,
    Cluster.info         = Cluster.info,
    Plot.cluster         = Plot.cluster,
    Cluster.merged       = Cluster.merged,
    Cluster.height       = Cluster.height,
    Species.cluster.phi  = Species.cluster.phi,
    species              = species,
    plots                = plots
  )
  class(res) <- c("cocktail", class(res))
  res
}
