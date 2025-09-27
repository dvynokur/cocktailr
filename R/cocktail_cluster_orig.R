#' Cocktail clustering (original version)
#'
#' Runs the original Bruelheide-like Cocktail clustering on a plots x species
#' matrix. Converts to presence/absence by default, computes phi-based merges,
#' and determines the per-cluster m using expected vs. observed frequencies.
#'
#' @inheritParams cocktail_cluster_new
#' @return A list (class "cocktail") with:
#' \itemize{
#'   \item Cluster.species  (clusters x species 0/1)
#'   \item Cluster.info     (k, m per cluster)
#'   \item Plot.cluster     (plots x clusters 0/1)
#'   \item Cluster.merged   (merge matrix like hclust)
#'   \item Cluster.height   (phi at merge)
#'   \item species, plots   (names)
#' }
#' @references Bruelheide (2000, 2016)
#' @export
cocktail_cluster_orig <- function(veg, to_pa = TRUE, progress = interactive()) {
  V <- as.matrix(veg)
  rownames(V) <- rownames(veg); colnames(V) <- colnames(veg)
  if (to_pa) V <- (V > 0) * 1L
  out <- .cocktail_core_orig(V, progress = progress)
  class(out) <- c("cocktail", class(out))
  out
}

# --- internal helpers ---------------------------------

.Expected.plot.freq_orig <- function(species_in_cluster, p_freq, q_freq) {
  k <- length(species_in_cluster)
  Exp <- array(0, k + 1)
  Exp_int <- array(0, k + 1)
  Exp_int[1] <- 1
  for (j in seq_len(k)) {
    s <- species_in_cluster[j]
    Exp[1] <- Exp_int[1] * q_freq[s]
    if (j > 1) {
      for (kk in 1:(j - 1)) {
        Exp[kk + 1] <- Exp_int[kk] * p_freq[s] + Exp_int[kk + 1] * q_freq[s]
      }
    }
    Exp[j + 1] <- Exp_int[j] * p_freq[s]
    for (kk in 1:(j + 1)) Exp_int[kk] <- Exp[kk]
  }
  Exp
}

.Compare.obs.exp.freq_orig <- function(Obs.freq, Exp.freq) {
  k <- length(Obs.freq) - 1L
  Cum.obs <- array(0, k + 1L)
  Cum.exp <- array(0, k + 1L)
  Cum.obs[k + 1L] <- Obs.freq[k + 1L]
  Cum.exp[k + 1L] <- Exp.freq[k + 1L]
  m <- 1L
  m.found <- -1L
  if (Cum.obs[k + 1L] > Cum.exp[k + 1L]) { m <- k; m.found <- 0L }
  for (j in seq(k, 1L, -1L)) {
    Cum.obs[j] <- Cum.obs[j + 1L] + Obs.freq[j]
    if (m.found == -1L && Cum.obs[j] > 0) m.found <- 0L
    Cum.exp[j] <- Cum.exp[j + 1L] + Exp.freq[j]
    if (j > 1L && m.found == 0L && Cum.exp[j] > Cum.obs[j]) { m <- j; m.found <- 1L }
  }
  m
}

# --- internal: clustering core -----------------------------

.cocktail_core_orig <- function(V_pa, progress = TRUE) {
  N <- nrow(V_pa); n <- ncol(V_pa)
  species <- colnames(V_pa); plots <- rownames(V_pa)

  p.freq <- colSums(V_pa) / N
  q.freq <- 1 - p.freq

  Cluster.species <- array(0L, c(n - 1L, n), dimnames = list(NULL, species))
  Cluster.info    <- array(0L, c(n - 1L, 2L), dimnames = list(seq_len(n - 1L), c("k","m")))
  Plot.cluster    <- array(0L, c(N, n - 1L), dimnames = list(plots, seq_len(n - 1L)))
  Cluster.merged  <- array(0L, c(n - 1L, 2L))
  Cluster.height  <- array(0,  n - 1L)

  veg2 <- V_pa
  n2 <- n
  i <- 0L
  name.last.cluster <- NULL

  pb <- if (progress) utils::txtProgressBar(min = 0, max = n - 1L, style = 3) else NULL
  on.exit(if (!is.null(pb)) close(pb), add = TRUE)

  while (i <= (n - 2L)) {
    phi.index <- vegan::designdist(
      t(veg2),
      method="(a*d-b*c)/sqrt((a+c)*(b+d)*(a+b)*(c+d))",
      terms="binary", alphagamma=FALSE, abcd=TRUE, "phi"
    )
    phi.index[is.na(phi.index)] <- 0

    allmax <- which(phi.index == max(phi.index))
    elements.per.col <- vapply(1:(n2 - 1L), function(x) (x - 1L) * n2 - ((x - 1L) * x / 2), numeric(1))
    multiple.max <- length(allmax)
    e1 <- integer(multiple.max)
    e2 <- integer(multiple.max)
    for (j in 1:multiple.max) {
      e1[j] <- max(which(elements.per.col < allmax[j]))
      e2[j] <- allmax[j] - elements.per.col[e1[j]] + e1[j]
    }

    compare1 <- as.vector(t(as.matrix(cbind(e1, e2))))
    if (anyDuplicated(compare1) > 2) {
      keep <- rep(TRUE, multiple.max)
      for (j in 2:multiple.max) {
        compare2 <- compare1[1:(2 * (j - 1))]
        k1 <- length(stats::na.omit(match(compare2, e1[j])))
        k2 <- length(stats::na.omit(match(compare2, e2[j])))
        if (k1 > 0 & k2 > 0) keep[j] <- FALSE
      }
      e1 <- e1[keep]
      e2 <- e2[keep]
      multiple.max <- length(e1)
    }

    i1 <- i + 1L
    for (j in 1:multiple.max) {
      i <- i + 1L
      Cluster.height[i] <- phi.index[allmax][j]

      if (substr(colnames(veg2)[e1[j]], 1, 2) == "c_") {
        cl1 <- as.integer(strsplit(colnames(veg2)[e1[j]], "_")[[1]][2])
        Cluster.merged[i, 1] <- cl1
        Cluster.species[i, which(Cluster.species[cl1, ] == 1)] <- 1L
      } else {
        f1 <- which(colnames(V_pa) == colnames(veg2)[e1[j]])
        Cluster.merged[i, 1] <- -f1
        Cluster.species[i, f1] <- 1L
      }

      if (substr(colnames(veg2)[e2[j]], 1, 2) == "c_") {
        cl2 <- as.integer(strsplit(colnames(veg2)[e2[j]], "_")[[1]][2])
        Cluster.merged[i, 2] <- cl2
        Cluster.species[i, which(Cluster.species[cl2, ] == 1)] <- 1L
      } else {
        f2 <- which(colnames(V_pa) == colnames(veg2)[e2[j]])
        Cluster.merged[i, 2] <- -f2
        Cluster.species[i, f2] <- 1L
      }

      # rename merged columns to new cluster id
      cn <- colnames(veg2)
      cn[cn == colnames(veg2)[e1[j]]] <- paste0("c_", i)
      cn[cn == colnames(veg2)[e2[j]]] <- paste0("c_", i)
      colnames(veg2) <- cn

      # k, m and plot assignment
      Cluster.info[i, "k"] <- sum(Cluster.species[i, ])
      spp_in_plot <- rowSums(V_pa[, which(Cluster.species[i, ] > 0), drop = FALSE])
      Obs <- matrix(0, Cluster.info[i, "k"] + 1L, dimnames = list(0:Cluster.info[i, "k"]))
      tt <- table(spp_in_plot)
      Obs[, 1] <- as.numeric(tt[match(rownames(Obs), names(tt))])
      Obs[is.na(Obs)] <- 0
      Exp <- .Expected.plot.freq_orig(which(Cluster.species[i, ] > 0), p.freq, q.freq) * N
      Cluster.info[i, "m"] <- .Compare.obs.exp.freq_orig(Obs, Exp)
      Plot.cluster[spp_in_plot >= Cluster.info[i, "m"], i] <- 1L
    }

    if (n2 == 3L) name.last.cluster <- colnames(veg2)[-c(e1[1], e2[1])]
    index.e1 <- !duplicated(colnames(veg2)[e1], fromLast = TRUE)
    index.e2 <- !duplicated(colnames(veg2)[e2], fromLast = TRUE)
    index.e <- seq(i1, i)[index.e1 & index.e2]
    veg2 <- cbind(veg2[, -unique(c(e1, e2)), drop = FALSE],
                  Plot.cluster[, index.e, drop = FALSE])
    n2 <- ncol(veg2)
    for (jj in seq_along(index.e)) {
      colnames(veg2)[n2 - length(index.e) + jj] <- paste0("c_", index.e[jj])
    }
    if (n2 == 2L && length(name.last.cluster) == 1L) {
      colnames(veg2)[1] <- name.last.cluster
    }

    if (sum(Plot.cluster[, i]) == N) {
      for (j in (i + 1L):(n - 1L)) {
        g1 <- ncol(veg2)
        cl1 <- as.integer(strsplit(colnames(veg2)[g1], "_")[[1]][2])
        Cluster.merged[j, 1] <- cl1
        g2 <- 1L
        cl2 <- as.integer(strsplit(colnames(veg2)[g2], "_")[[1]][2])
        Cluster.merged[j, 2] <- cl2
        Cluster.height[j] <- 0
        Plot.cluster[, j] <- 1L
        Cluster.info[j, "k"] <- sum(Cluster.species[j, ])
        Cluster.info[j, "m"] <- 1L
        Cluster.species[j, which(Cluster.species[cl1, ] == 1)] <- 1L
        Cluster.species[j, which(Cluster.species[cl2, ] == 1)] <- 1L
        veg2 <- cbind(veg2, Plot.cluster[, j, drop = FALSE])
        colnames(veg2)[ncol(veg2)] <- paste0("c_", j)
        veg2 <- veg2[, -c(g1, g2), drop = FALSE]
      }
      i <- j
    }

    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    if (i >= (n - 1L)) break
  }

  list(
    Cluster.species = Cluster.species,
    Cluster.info    = Cluster.info,
    Plot.cluster    = Plot.cluster,
    Cluster.merged  = Cluster.merged,
    Cluster.height  = Cluster.height,
    species = species, plots = plots
  )
}
