#' Cocktail clustering (sparse matrix)
#'
#' Fast Cocktail agglomeration for a **plots × species** table.
#'
#' @description
#' This implementation:
#' - **Binarizes** the input: values > 0 become 1; values ≤ 0 or `NA` become 0.
#' - **Drops empty species** (all-zero columns) before clustering.
#' - Computes the association coefficient (“phi”) each round from one sparse
#'   crossproduct for speed and exactness.
#' - Uses a **fixed, reproducible tie order**: when several pairs share the same
#'   maximum phi at a step, they are processed in the same order that R fills the
#'   lower-triangular distance matrix (scan by increasing column, then row).
#'
#' @param vegmatrix A matrix or data frame with **plots in rows** and **species in columns**.
#' @param progress  Logical; show a text progress bar (default `TRUE`).
#'
#' @return
#' A list of class `"cocktail"` with:
#' \itemize{
#'   \item `Cluster.species` — integer matrix (n_merges × n_species): species membership per merge.
#'   \item `Cluster.merged`  — integer matrix (n_merges × 2): left/right children per merge
#'         (negative = original species index; positive = earlier merge index).
#'   \item `Cluster.height`  — numeric vector length n_merges: phi at each merge.
#'   \item `Plot.cluster`    — integer matrix (n_plots × n_merges): plot membership per merge.
#'   \item `species`         — character vector of species names kept after cleaning.
#'   \item `plots`           — character vector of plot names.
#' }
#'
#' @details
#' Binarization and removal of empty species happen internally and only affect the
#' set of columns that contribute to clustering. All returned components are aligned
#' to the species that had at least one presence after cleaning.
#'
#' @seealso \code{\link{plot_cocktail}}, \code{\link{cocktail_fuzzy}}, \code{\link{assign_releves}}
#'
#' @examples
#' vm <- matrix(c(1,0,1,
#'                0,1,0,
#'                1,1,0),
#'              nrow = 3, byrow = TRUE,
#'              dimnames = list(paste0("plot", 1:3),
#'                              c("sp1","sp2","sp3")))
#' res <- cocktail_cluster(vm, progress = FALSE)
#' names(res)
#'
#' @import Matrix
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export

cocktail_cluster <- function(
    vegmatrix,
    progress = TRUE
) {
  ## ---- input checks & setup ----
  if (!is.matrix(vegmatrix) && !is.data.frame(vegmatrix)) {
    stop("vegmatrix must be a matrix or data.frame with plots in rows and species in columns.")
  }
  vm <- as.matrix(vegmatrix)

  vm[is.na(vm)] <- 0
  vm <- (vm > 0) * 1L

  plots   <- rownames(vm); if (is.null(plots))   plots   <- as.character(seq_len(nrow(vm)))
  species <- colnames(vm); if (is.null(species)) species <- as.character(seq_len(ncol(vm)))

  # Drop empty species columns
  keep <- colSums(vm) > 0L
  dropped_species <- character(0)
  if (!all(keep)) {
    dropped_species <- species[!keep]
    vm      <- vm[, keep, drop = FALSE]
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

  # idx(i,j; m) = (i-1)*m - ((i-1)*i)/2 + (j - i),  for 1<=i<j<=m
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
  Plot.cluster    <- matrix(0L, N, n - 1L, dimnames = list(plots, NULL))
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
    if (!nrow(maxres$pairs) || !(is.finite(maxres$max_phi) && maxres$max_phi > 0)) {
      maxres <- phi_max_pairs_fallback_distorder_(X)
      if (!nrow(maxres$pairs)) break
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

      Plot.cluster[species_in_plot >= Cluster.info[i, 2L], i] <- 1L
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
        Matrix::Matrix(Plot.cluster[, kid] > 0L, sparse = TRUE)
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
    if (sum(Plot.cluster[, i]) == N) {
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
          Plot.cluster[, j] <- 1L
          Cluster.info[j, 1L] <- sum(Cluster.species[j, ])
          Cluster.info[j, 2L] <- 1L

          Cluster.species[j, Cluster.species[cl1, ] == 1L] <- 1L
          Cluster.species[j, Cluster.species[cl2, ] == 1L] <- 1L

          X <- cbind(X, Matrix::Matrix(Plot.cluster[, j] > 0L, sparse = TRUE))
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

  res <- list(
    Cluster.species = Cluster.species,
    Cluster.info    = Cluster.info,
    Plot.cluster    = Plot.cluster,
    Cluster.merged  = Cluster.merged,
    Cluster.height  = Cluster.height,
    species         = species,
    plots           = plots#,
    # dropped_species = dropped_species
  )
  class(res) <- c("cocktail", class(res))
  res
}
