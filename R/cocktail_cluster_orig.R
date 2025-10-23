#' Original Cocktail clustering (reference implementation)
#'
#' Run the published Cocktail agglomeration algorithm on a plots × species table.
#'
#' @description Temporary reference implementation to compare with `cocktail_cluster()`.
#'
#' @param vegmatrix Matrix/data.frame, plots in rows, species in columns.
#' @param progress Logical; show a text progress bar.
#' @return A list with Cluster.species, Cluster.info, Plot.cluster, Cluster.merged, Cluster.height, species, plots.
#' @examples
#' \donttest{
#' m <- as.matrix(vegan::dune); m(m > 0) <- 1L
#' res <- cocktail_cluster_orig(m, progress = FALSE)
#' }
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom vegan designdist

cocktail_cluster_orig <- function(vegmatrix, progress = TRUE) {
  # ---- input & cleaning -----------------------------------------------------
  if (!is.matrix(vegmatrix) && !is.data.frame(vegmatrix)) {
    stop("vegmatrix must be a matrix or data.frame with plots in rows and species in columns.")
  }
  vegmatrix <- as.matrix(vegmatrix)
  vegmatrix[is.na(vegmatrix)] <- 0
  vegmatrix <- (vegmatrix > 0) * 1L

  plots   <- rownames(vegmatrix); if (is.null(plots))   plots   <- as.character(seq_len(nrow(vegmatrix)))
  species <- colnames(vegmatrix); if (is.null(species)) species <- as.character(seq_len(ncol(vegmatrix)))

  # align with cocktail_cluster(): drop empty species
  keep <- colSums(vegmatrix) > 0L
  if (!all(keep)) {
    vegmatrix <- vegmatrix[, keep, drop = FALSE]
    species   <- species[keep]
  }

  N <- nrow(vegmatrix)
  n <- ncol(vegmatrix)
  if (n < 2L) stop("After dropping empty species, need at least 2 species (columns).")
  if (N < 1L) stop("Need at least 1 plot (row).")

  # global freqs
  col.sums <- colSums(vegmatrix)
  row.sums <- rowSums(vegmatrix)
  p.freq <- col.sums / N
  q.freq <- 1 - p.freq

  # ---- helpers (verbatim logic) --------------------------------------------
  Expected.plot.freq <- function(species.in.cluster) {
    No.of.spec.in.cluster <- length(species.in.cluster)
    Exp.no.of.plots <- array(0, No.of.spec.in.cluster + 1L,
                             dimnames = list(seq(0, No.of.spec.in.cluster, 1)))
    Exp.no.of.plots.inter <- array(0, No.of.spec.in.cluster + 1L,
                                   dimnames = list(seq(0, No.of.spec.in.cluster, 1)))
    Exp.no.of.plots.inter[1L] <- 1
    for (j in 1L:No.of.spec.in.cluster) {
      Exp.no.of.plots[1L] <- Exp.no.of.plots.inter[1L] * q.freq[species.in.cluster[j]]
      if (j > 1L) {
        for (k in 1L:(j - 1L)) {
          Exp.no.of.plots[k + 1L] <- Exp.no.of.plots.inter[k]     * p.freq[species.in.cluster[j]] +
            Exp.no.of.plots.inter[k + 1L] * q.freq[species.in.cluster[j]]
        }
      }
      Exp.no.of.plots[j + 1L] <- Exp.no.of.plots.inter[j] * p.freq[species.in.cluster[j]]
      for (k in 1L:(j + 1L)) {
        Exp.no.of.plots.inter[k] <- Exp.no.of.plots[k]
      }
    }
    Exp.no.of.plots
  }

  Compare.obs.exp.freq <- function(Obs.freq, Exp.freq) {
    No.of.spec.in.cluster <- length(Obs.freq) - 1L
    Cum.obs.no.of.plots <- array(0, No.of.spec.in.cluster + 1L,
                                 dimnames = list(seq(0, No.of.spec.in.cluster, 1)))
    Cum.exp.no.of.plots <- array(0, No.of.spec.in.cluster + 1L,
                                 dimnames = list(seq(0, No.of.spec.in.cluster, 1)))
    Cum.obs.no.of.plots[No.of.spec.in.cluster + 1L] <- Obs.freq[No.of.spec.in.cluster + 1L]
    Cum.exp.no.of.plots[No.of.spec.in.cluster + 1L] <- Exp.freq[No.of.spec.in.cluster + 1L]
    m <- 1L
    m.found <- -1L
    if (Cum.obs.no.of.plots[No.of.spec.in.cluster + 1L] >
        Cum.exp.no.of.plots[No.of.spec.in.cluster + 1L]) {
      m <- No.of.spec.in.cluster
      m.found <- 0L
    }
    for (j in seq(No.of.spec.in.cluster, 1L, -1L)) {
      Cum.obs.no.of.plots[j] <- Cum.obs.no.of.plots[j + 1L] + Obs.freq[j]
      if (m.found == -1L && Cum.obs.no.of.plots[j] > 0) m.found <- 0L
      Cum.exp.no.of.plots[j] <- Cum.exp.no.of.plots[j + 1L] + Exp.freq[j]
      if (j > 1L && m.found == 0L && Cum.exp.no.of.plots[j] > Cum.obs.no.of.plots[j]) {
        m <- j
        m.found <- 1L
      }
    }
    m
  }

  # ---- outputs --------------------------------------------------------------
  Cluster.species <- array(0L, c(n - 1L, n))
  dimnames(Cluster.species)[[2]] <- species
  Cluster.info   <- array(0L, c(n - 1L, 2L), dimnames = list(c(1:(n - 1L)), c("k", "m")))
  Plot.cluster   <- array(0L, c(N, n - 1L), dimnames = list(plots, c(1:(n - 1L))))
  Cluster.merged <- array(0L, c(n - 1L, 2L))
  Cluster.height <- array(0, n - 1L)

  vegmatrix2 <- vegmatrix
  n2 <- n
  i <- 0L
  multiple.max <- 1L
  name.last.cluster <- NULL

  pb <- if (isTRUE(progress)) utils::txtProgressBar(min = 0, max = n - 1L, style = 3) else NULL
  on.exit({ if (!is.null(pb)) close(pb) }, add = TRUE)

  # ---- clustering loop ------------------------------------------------------
  while (i <= (n - 2L)) {
    # phi distance (actually similarity as in the original code)
    phi.index <- vegan::designdist(
      t(vegmatrix2),
      method = "(a*d-b*c)/sqrt((a+c)*(b+d)*(a+b)*(c+d))",
      terms = c("binary"),
      alphagamma = FALSE,
      abcd = TRUE
    )
    phi.index[is.na(phi.index)] <- 0

    allmax <- which(phi.index == max(phi.index))
    elements.per.col <- array(0, n2 - 1L)
    for (x in 1:(n2 - 1L)) {
      elements.per.col[x] <- (x - 1L) * n2 - ((x - 1L) * x / 2)
    }
    multiple.max <- length(allmax)
    e1 <- array(0L, multiple.max)
    e2 <- array(0L, multiple.max)
    for (j in 1:multiple.max) {
      e1[j] <- max(which(elements.per.col < allmax[j]))
      e2[j] <- allmax[j] - elements.per.col[e1[j]] + e1[j]
    }

    # circularity filter
    compare1 <- as.vector(t(as.matrix(cbind(e1, e2))))
    if (anyDuplicated(compare1) > 2L) {
      tobekept <- rep(TRUE, multiple.max)
      for (j in 2:multiple.max) {
        compare2 <- compare1[1:(2 * (j - 1))]
        k1 <- length(match(compare2, e1[j])[!is.na(match(compare2, e1[j]))])
        k2 <- length(match(compare2, e2[j])[!is.na(match(compare2, e2[j]))])
        if (k1 > 0 && k2 > 0) tobekept[j] <- FALSE
      }
      e1 <- e1[tobekept]
      e2 <- e2[tobekept]
      multiple.max <- length(e1)
      if (multiple.max == 0L) next
    }

    i1 <- i + 1L
    for (j in 1:multiple.max) {
      i <- i + 1L
      Cluster.height[i] <- phi.index[allmax][j]

      if (substr(colnames(vegmatrix2)[e1[j]], 1, 2) == "c_") {
        cluster.no.1 <- as.numeric(strsplit(colnames(vegmatrix2)[e1[j]], split = "_")[[1]][2])
        Cluster.merged[i, 1] <- cluster.no.1
        Cluster.species[i, which(Cluster.species[cluster.no.1, ] == 1L)] <- 1L
      } else {
        f1 <- which(species == colnames(vegmatrix2)[e1[j]])
        Cluster.merged[i, 1] <- -f1
        Cluster.species[i, f1] <- 1L
      }
      if (substr(colnames(vegmatrix2)[e2[j]], 1, 2) == "c_") {
        cluster.no.2 <- as.numeric(strsplit(colnames(vegmatrix2)[e2[j]], split = "_")[[1]][2])
        Cluster.merged[i, 2] <- cluster.no.2
        Cluster.species[i, which(Cluster.species[cluster.no.2, ] == 1L)] <- 1L
      } else {
        f2 <- which(species == colnames(vegmatrix2)[e2[j]])
        Cluster.merged[i, 2] <- -f2
        Cluster.species[i, f2] <- 1L
      }

      # rename endpoints of current pair to new cluster name
      colnames(vegmatrix2)[c(e1, e2)][colnames(vegmatrix2)[c(e1, e2)] == colnames(vegmatrix2)[e1[j]]] <-
        paste0("c_", i)
      colnames(vegmatrix2)[c(e1, e2)][colnames(vegmatrix2)[c(e1, e2)] == colnames(vegmatrix2)[e2[j]]] <-
        paste0("c_", i)

      # k, observed/expected plot freqs, m, plot assignment
      Cluster.info[i, 1] <- sum(Cluster.species[i, ])
      spp_idx <- which(Cluster.species[i, ] > 0L)
      species.in.plot <- rowSums(vegmatrix[, spp_idx, drop = FALSE])
      species.in.plot.table <- as.matrix(table(species.in.plot))
      Obs.plot.freq <- matrix(0, Cluster.info[i, 1] + 1L,
                              dimnames = list(0:Cluster.info[i, 1]))
      index <- match(dimnames(Obs.plot.freq)[[1]],
                     dimnames(table(species.in.plot))$species.in.plot)
      Obs.plot.freq[, 1] <- as.numeric(table(species.in.plot))[index]
      Obs.plot.freq <- ifelse(is.na(Obs.plot.freq), 0, Obs.plot.freq)
      Exp.plot.freq <- Expected.plot.freq(spp_idx) * N
      Cluster.info[i, 2] <- Compare.obs.exp.freq(Obs.plot.freq, Exp.plot.freq)
      Plot.cluster[species.in.plot >= Cluster.info[i, 2], i] <- 1L
    }

    if (n2 == 3L) {
      name.last.cluster <- colnames(vegmatrix2)[-c(e1[1], e2[1])]
    }

    index.e1 <- !duplicated(colnames(vegmatrix2)[e1], fromLast = TRUE)
    index.e2 <- !duplicated(colnames(vegmatrix2)[e2], fromLast = TRUE)
    index.e <- c(i1:i)[index.e1 & index.e2]

    vegmatrix2 <- cbind(vegmatrix2[, -unique(c(e1, e2)), drop = FALSE],
                        Plot.cluster[, index.e, drop = FALSE])
    n2 <- ncol(vegmatrix2)
    for (j in seq_along(index.e)) {
      colnames(vegmatrix2)[n2 - length(index.e) + j] <- paste0("c_", index.e[j])
    }
    if (n2 == 2L) {
      colnames(vegmatrix2)[1] <- name.last.cluster
    }

    # stop criterion: all plots merged → fill remaining merges at phi = 0
    if (sum(Plot.cluster[, i]) == N) {
      for (j in (i + 1):(n - 1L)) {
        g1 <- ncol(vegmatrix2)
        cluster.no.1 <- as.numeric(strsplit(colnames(vegmatrix2)[g1], split = "_")[[1]][2])
        Cluster.merged[j, 1] <- cluster.no.1
        g2 <- 1L
        cluster.no.2 <- as.numeric(strsplit(colnames(vegmatrix2)[g2], split = "_")[[1]][2])
        Cluster.merged[j, 2] <- cluster.no.2
        Cluster.height[j] <- 0
        Plot.cluster[, j] <- 1L
        Cluster.info[j, 1] <- sum(Cluster.species[j, ])
        Cluster.info[j, 2] <- 1L
        Cluster.species[j, which(Cluster.species[cluster.no.1, ] == 1L)] <- 1L
        Cluster.species[j, which(Cluster.species[cluster.no.2, ] == 1L)] <- 1L
        vegmatrix2 <- cbind(vegmatrix2, Plot.cluster[, j, drop = FALSE])
        colnames(vegmatrix2)[ncol(vegmatrix2)] <- paste0("c_", j)
        vegmatrix2 <- vegmatrix2[, -c(g1, g2), drop = FALSE]
        n2 <- ncol(vegmatrix2)
      }
      i <- j
    }

    if (!is.null(pb)) utils::setTxtProgressBar(pb, min(i, n - 1L))
  }

  list(
    Cluster.species = Cluster.species,
    Cluster.info    = Cluster.info,
    Plot.cluster    = Plot.cluster,
    Cluster.merged  = Cluster.merged,
    Cluster.height  = Cluster.height,
    species         = species,
    plots           = plots
  )
}
