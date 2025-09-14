#' Cocktail clustering (sparse, parity with original, bounded merges)
#'
#' Runs the Cocktail hierarchical agglomeration on a plots x species matrix.
#' Each round uses a single sparse crossproduct to compute co-occurrences and
#' φ from (a,b,c,d) via the same (ad−bc)/sqrt(...) form as vegan::designdist.
#' Ties are ordered to match the lower-triangle "dist" indexing. When the
#' current max φ is nonpositive, an optional fallback evaluates the round
#' with vegan::designdist for exact parity (including a=0 pairs).
#'
#' @param vegmatrix Numeric matrix or data.frame; plots in rows, species in columns.
#' @param binarize Logical; convert to presence/absence (default TRUE).
#' @param progress Logical; show a text progress bar (default TRUE).
#' @param fallback_when_nonpos Logical; if TRUE, when max φ ≤ 0 compute that round
#'   with vegan::designdist for exact parity. Default TRUE.
#' @param fallback_max_cols Integer; only use the fallback if the current number
#'   of columns is ≤ this value (prevents stalls on huge rounds). Default 600.
#'
#' @return A list of class \code{"cocktail"} with:
#' \itemize{
#'   \item \code{Cluster.species}: (n-1) x n integer matrix; species membership per cluster
#'   \item \code{Cluster.info}:   (n-1) x 2 integer matrix; columns \code{k} (size), \code{m} (min spp per plot)
#'   \item \code{Plot.cluster}:    N x (n-1) integer matrix; plot assignment per cluster
#'   \item \code{Cluster.merged}:  (n-1) x 2 integer matrix; negative = species id, positive = cluster id
#'   \item \code{Cluster.height}:  numeric vector of length n-1; phi at which each merge formed
#'   \item \code{species}, \code{plots}: character vectors with names from input
#' }
#'
#' @import Matrix
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom vegan designdist
#' @export
cocktail_cluster_new <- function(
    vegmatrix,
    binarize = TRUE,
    progress = TRUE,
    fallback_when_nonpos = TRUE,
    fallback_max_cols = 600L
) {
  ## ---- input checks and setup ----
  if (!is.matrix(vegmatrix) && !is.data.frame(vegmatrix)) {
    stop("vegmatrix must be a matrix or data.frame with plots in rows and species in columns.")
  }
  vm <- as.matrix(vegmatrix)
  storage.mode(vm) <- "integer"
  if (isTRUE(binarize)) vm[] <- as.integer(vm > 0L)

  plots   <- rownames(vm); if (is.null(plots))   plots   <- as.character(seq_len(nrow(vm)))
  species <- colnames(vm); if (is.null(species)) species <- as.character(seq_len(ncol(vm)))

  N <- nrow(vm); n <- ncol(vm)
  if (n < 2L) stop("Need at least 2 species (columns).")
  if (N < 1L) stop("Need at least 1 plot (row).")

  # sparse working matrix (plots x species/clusters)
  X0 <- Matrix::Matrix(vm, sparse = TRUE)  # dgCMatrix
  rm(vm)

  # global frequencies for Expected.plot.freq (original recursion)
  p.freq <- as.numeric(Matrix::colSums(X0)) / N
  q.freq <- 1 - p.freq

  ## ---- helpers ----
  Expected.plot.freq_ <- function(species.in.cluster) {
    No <- length(species.in.cluster)
    Exp.no       <- array(0, No + 1L)
    Exp.no.inter <- array(0, No + 1L)
    Exp.no.inter[1L] <- 1
    for (j in 1L:No) {
      s <- species.in.cluster[j]
      Exp.no[1L] <- Exp.no.inter[1L] * q.freq[s]
      if (j > 1L) {
        for (k in 1L:(j - 1L)) {
          Exp.no[k + 1L] <- Exp.no.inter[k]     * p.freq[s] +
            Exp.no.inter[k + 1L] * q.freq[s]
        }
      }
      Exp.no[j + 1L] <- Exp.no.inter[j] * p.freq[s]
      for (k in 1L:(j + 1L)) Exp.no.inter[k] <- Exp.no[k]
    }
    Exp.no
  }

  Compare.obs.exp.freq_ <- function(Obs.freq, Exp.freq) {
    Obs <- if (is.matrix(Obs.freq)) as.vector(Obs.freq[, 1L]) else as.vector(Obs.freq)
    No <- length(Obs) - 1L
    Cum.obs <- array(0, No + 1L)
    Cum.exp <- array(0, No + 1L)
    Cum.obs[No + 1L] <- Obs[No + 1L]
    Cum.exp[No + 1L] <- Exp.freq[No + 1L]
    m <- 1L
    m.found <- -1L
    if (Cum.obs[No + 1L] > Cum.exp[No + 1L]) { m <- No; m.found <- 0L }
    for (j in No:1L) {
      Cum.obs[j] <- Cum.obs[j + 1L] + Obs[j]
      if (m.found == -1L && Cum.obs[j] > 0) m.found <- 0L
      Cum.exp[j] <- Cum.exp[j + 1L] + Exp.freq[j]
      if (j > 1L && m.found == 0L && Cum.exp[j] > Cum.obs[j]) { m <- j; m.found <- 1L }
    }
    m
  }

  # lower-triangle "dist" index to mimic vegan's tie ordering (1-based i<j)
  .dist_index_ <- function(i, j, m) (i - 1L) * m - ((i - 1L) * i) %/% 2L + (j - i)

  # φ by sparse crossproduct (ad−bc form), with ties ordered like dist
  phi_max_pairs_crossprod_ <- function(X) {
    stopifnot(inherits(X, "dgCMatrix"))
    N  <- nrow(X); m <- ncol(X)
    if (m < 2L) return(list(max_phi = 0, pairs = matrix(integer(0), ncol = 2)))

    p <- as.numeric(Matrix::colSums(X))

    # all co-occurrences at once (symmetric; keep upper triangle only)
    A <- Matrix::crossprod(X)      # dsCMatrix (upper triangle stored)
    A <- Matrix::triu(A, k = 1L)   # drop diagonal
    A <- as(A, "dgCMatrix")
    if (length(A@x) == 0L) {
      if (m >= 2L) return(list(max_phi = 0, pairs = cbind(e1 = 1L, e2 = 2L)))
      return(list(max_phi = 0, pairs = matrix(integer(0), ncol = 2)))
    }

    # indices aligned with A@x (i<j by construction)
    j_idx <- rep.int(seq_len(m), diff(A@p))
    i_idx <- A@i + 1L
    a     <- A@x

    # exact counts as in vegan::designdist(abcd=TRUE)
    b <- p[i_idx] - a
    c <- p[j_idx] - a
    d <- N - a - b - c

    num <- a * d - b * c
    den <- sqrt((a + b) * (c + d) * (a + c) * (b + d))
    den[den == 0] <- NA_real_

    phi <- num / den
    phi[!is.finite(phi)] <- 0

    max_phi <- max(phi, na.rm = TRUE)
    if (!is.finite(max_phi)) max_phi <- 0

    keep <- which(phi == max_phi)  # match original's exact equality
    e1 <- i_idx[keep]; e2 <- j_idx[keep]
    ord <- order(.dist_index_(e1, e2, m))
    list(max_phi = max_phi, pairs = cbind(e1 = e1[ord], e2 = e2[ord]))
  }

  # One-round φ finder with optional fallback when max φ ≤ 0
  phi_max_pairs_round_ <- function(X) {
    fast <- phi_max_pairs_crossprod_(X)
    if (!fallback_when_nonpos) return(fast)
    if (fast$max_phi > 0) return(fast)

    m <- ncol(X)
    if (m > fallback_max_cols) return(fast)  # avoid stalls on very large rounds

    # exact parity fallback (includes a=0 pairs; NA→0)
    M <- as.matrix(X)
    phi_full <- vegan::designdist(
      t(M),
      method="(a*d-b*c)/sqrt((a+c)*(b+d)*(a+b)*(c+d))",
      terms="binary", alphagamma=FALSE, abcd=TRUE, "phi"
    )
    PM <- as.matrix(phi_full); PM[is.na(PM)] <- 0
    maxv <- max(PM)
    e1 <- integer(0); e2 <- integer(0)
    for (ii in seq_len(m - 1L)) {
      jj <- (ii + 1L):m
      hits <- jj[PM[ii, jj] == maxv]
      if (length(hits)) { e1 <- c(e1, rep.int(ii, length(hits))); e2 <- c(e2, hits) }
    }
    list(max_phi = maxv, pairs = cbind(e1 = e1, e2 = e2))
  }

  ## ---- outputs ----
  Cluster.species <- matrix(0L, n - 1L, n, dimnames = list(NULL, species))
  Cluster.info    <- matrix(0L, n - 1L, 2L, dimnames = list(NULL, c("k","m")))
  Plot.cluster    <- matrix(0L, N, n - 1L, dimnames = list(plots, NULL))
  Cluster.merged  <- matrix(0L, n - 1L, 2L)
  Cluster.height  <- numeric(n - 1L)

  ## ---- working state ----
  X         <- X0
  col_names <- colnames(X0); if (is.null(col_names)) col_names <- species
  i <- 0L

  pb <- if (isTRUE(progress)) utils::txtProgressBar(min = 0, max = n - 1L, style = 3) else NULL
  on.exit({ if (!is.null(pb)) close(pb) }, add = TRUE)

  ## ---- agglomeration loop ----
  while (i <= (n - 2L)) {

    # find all pairs at current max φ
    maxres <- phi_max_pairs_round_(X)
    if (!nrow(maxres$pairs)) break
    e1 <- maxres$pairs[, "e1"]
    e2 <- maxres$pairs[, "e2"]
    multiple.max <- length(e1)

    # anti-circularization as in original
    compare1 <- as.vector(t(cbind(e1, e2)))
    if (anyDuplicated(compare1) > 2L) {
      tobekept <- rep(TRUE, multiple.max)
      for (jj in 2L:multiple.max) {
        compare2 <- compare1[1L:(2L * (jj - 1L))]
        k1 <- sum(!is.na(match(compare2, e1[jj])))
        k2 <- sum(!is.na(match(compare2, e2[jj])))
        if (k1 > 0L & k2 > 0L) tobekept[jj] <- FALSE
      }
      e1 <- e1[tobekept]; e2 <- e2[tobekept]
      multiple.max <- length(e1)
      if (multiple.max == 0L) next
    }

    # cap to avoid writing past (n-1) rows
    remaining_merges <- (n - 1L) - i
    if (remaining_merges <= 0L) break
    if (multiple.max > remaining_merges) {
      e1 <- e1[seq_len(remaining_merges)]
      e2 <- e2[seq_len(remaining_merges)]
      multiple.max <- remaining_merges
    }

    i1 <- i + 1L

    # fuse all maxima at this height
    for (jj in seq_len(multiple.max)) {
      if (i >= (n - 1L)) break
      i <- i + 1L
      Cluster.height[i] <- maxres$max_phi

      # left element
      left_name <- col_names[e1[jj]]
      if (startsWith(left_name, "c_")) {
        cl1 <- as.integer(sub("c_", "", left_name))
        Cluster.merged[i, 1L] <- cl1
        Cluster.species[i, Cluster.species[cl1, ] == 1L] <- 1L
      } else {
        f1 <- match(left_name, species)
        Cluster.merged[i, 1L] <- -f1
        Cluster.species[i, f1] <- 1L
      }

      # right element
      right_name <- col_names[e2[jj]]
      if (startsWith(right_name, "c_")) {
        cl2 <- as.integer(sub("c_", "", right_name))
        Cluster.merged[i, 2L] <- cl2
        Cluster.species[i, Cluster.species[cl2, ] == 1L] <- 1L
      } else {
        f2 <- match(right_name, species)
        Cluster.merged[i, 2L] <- -f2
        Cluster.species[i, f2] <- 1L
      }

      # assign new cluster name to both fused columns (NAME-BASED, across all columns)
      newc <- paste0("c_", i)
      old1 <- col_names[e1[jj]]
      old2 <- col_names[e2[jj]]
      col_names[col_names == old1] <- newc
      col_names[col_names == old2] <- newc

      # k
      k <- sum(Cluster.species[i, ])
      Cluster.info[i, 1L] <- k

      # observed plot-frequency (fast)
      spp_idx <- which(Cluster.species[i, ] > 0L)
      species_in_plot <- as.integer(Matrix::rowSums(X0[, spp_idx, drop = FALSE]))
      Obs_plot_freq <- tabulate(species_in_plot + 1L, nbins = k + 1L)

      # expected plot-frequency (×N)
      Exp_plot_freq <- Expected.plot.freq_(spp_idx) * N

      # m
      Cluster.info[i, 2L] <- Compare.obs.exp.freq_(Obs_plot_freq, Exp_plot_freq)

      # assign plots to cluster
      Plot.cluster[species_in_plot >= Cluster.info[i, 2L], i] <- 1L
    }

    # rebuild working matrix: drop fused cols and append the new cluster columns
    if (multiple.max > 0L) {
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
    }

    # fast tail: if the most recent cluster covers all plots, finish with φ=0 merges
    if (i > 0L && sum(Plot.cluster[, i]) == N) {
      for (jj in (i + 1L):(n - 1L)) {
        g1 <- ncol(X)
        cl1 <- as.integer(sub("c_", "", col_names[g1]))
        Cluster.merged[jj, 1L] <- cl1
        g2 <- 1L
        cl2 <- as.integer(sub("c_", "", col_names[g2]))
        Cluster.merged[jj, 2L] <- cl2
        Cluster.height[jj] <- 0
        Plot.cluster[, jj] <- 1L
        Cluster.info[jj, 1L] <- sum(Cluster.species[jj, ])
        Cluster.info[jj, 2L] <- 1L
        Cluster.species[jj, Cluster.species[cl1, ] == 1L] <- 1L
        Cluster.species[jj, Cluster.species[cl2, ] == 1L] <- 1L
        newc <- Matrix::Matrix(Plot.cluster[, jj] > 0L, sparse = TRUE)
        colnames(newc) <- paste0("c_", jj)
        X <- cbind(X, newc)
        X <- X[, -c(g1, g2), drop = FALSE]
        col_names <- c(col_names, colnames(newc))
        col_names <- col_names[-c(g1, g2)]
      }
      i <- n - 1L
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
    plots           = plots
  )
  class(res) <- c("cocktail", class(res))
  res
}
