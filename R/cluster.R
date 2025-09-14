#' Cocktail clustering
#'
#' Runs the Cocktail hierarchical agglomeration on a plots x species matrix,
#' returning the same set of objects as in the original code:
#' \itemize{
#'   \item \code{Cluster.species}: (n-1) x n binary matrix; species membership per cluster
#'   \item \code{Cluster.info}:   (n-1) x 2 matrix with columns \code{k} (size), \code{m} (min spp per plot)
#'   \item \code{Plot.cluster}:    N x (n-1) binary matrix; plot assignment per cluster
#'   \item \code{Cluster.merged}:  (n-1) x 2 matrix; negative = species id, positive = cluster id
#'   \item \code{Cluster.height}:  length n-1 numeric; phi value at which each merge formed
#'   \item \code{species}:         character vector of species (column names of the input)
#'   \item \code{plots}:           character vector of plot IDs (row names of the input)
#' }
#'
#' @param vegmatrix numeric matrix or data.frame, \strong{plots in rows, species in columns}.
#'   If you currently have species in rows, call \code{t()} before passing here.
#' @param binarize logical; if TRUE (default) convert positive abundances to presence/absence.
#' @param progress logical; show a text progress bar (default TRUE).
#' @param save_dir optional directory; if provided and \code{save_csv} or \code{save_rds} are TRUE,
#'   files will be written there (created if missing).
#' @param save_csv logical; write the four CSVs \code{Cluster.species.csv}, \code{Cluster.info.csv},
#'   \code{Plot.cluster.csv}, \code{Cluster.merged.csv} (default FALSE).
#' @param save_rds logical; write one RDS bundle with the full returned list (default FALSE).
#' @param rds_filename filename for the RDS bundle (default \code{"cocktail_cluster.rds"}).
#' @param rds_compress one of \code{"xz"}, \code{"gzip"}, \code{"bzip2"}, \code{"none"} (case-insensitive).
#'   If an invalid value is supplied, it silently falls back to \code{"none"} and warns.
#'
#' @details
#' The algorithm computes the pairwise \emph{phi} coefficient on species,
#' merges the most similar pairs (handling ties), and for each new group
#' determines the \code{m} threshold (minimum species per plot) by comparing
#' observed vs. expected plot-frequency distributions (Bruelheide 2000).
#'
#' Memory tips (large data): set a larger workspace and disable GUI graphics, e.g.
#' \preformatted{
#'   options(expressions = 5e5)
#'   options(device = "pdf")    # avoid GUI devices on Windows
#' }
#'
#' \strong{Forgot to assign the result?} If you call \code{cocktail_cluster(x)} without
#' assignment, you can still retrieve the last returned value via \code{.Last.value}
#' provided you haven't run another command that changed it.
#'
#' @return A list with elements described above; class \code{"cocktail"} is added.
#' @export
#' @importFrom vegan designdist
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv

cocktail_cluster <- function(
    vegmatrix,
    binarize   = TRUE,
    progress   = TRUE,
    save_dir   = NULL,
    save_csv   = FALSE,
    save_rds   = FALSE,
    rds_filename = "cocktail_cluster.rds",
    rds_compress = c("xz","gzip","bzip2","none")
) {
  ## -------- Input checks --------
  if (!is.matrix(vegmatrix) && !is.data.frame(vegmatrix)) {
    stop("vegmatrix must be a matrix or data.frame (plots in rows, species in columns).")
  }
  vm <- as.matrix(vegmatrix)
  storage.mode(vm) <- "numeric"

  plots   <- rownames(vm)
  species <- colnames(vm)
  if (is.null(plots))   plots   <- as.character(seq_len(nrow(vm)))
  if (is.null(species)) species <- as.character(seq_len(ncol(vm)))

  if (binarize) vm[vm > 0] <- 1

  N <- nrow(vm)  # plots
  n <- ncol(vm)  # species
  if (n < 2L) stop("Need at least 2 species (columns).")
  if (N < 1L) stop("Need at least 1 plot (row).")

  ## -------- Frequencies (global) --------
  col.sums <- colSums(vm)
  p.freq <- col.sums / N
  q.freq <- 1 - p.freq

  ## -------- Helper functions (as in original) --------
  Expected.plot.freq <- function(species.in.cluster) {
    No <- length(species.in.cluster)
    Exp.no      <- array(0, No + 1, dimnames = list(seq(0, No, 1)))
    Exp.no.inter<- array(0, No + 1, dimnames = list(seq(0, No, 1)))
    Exp.no.inter[1] <- 1
    for (j in 1:No) {
      s <- species.in.cluster[j]
      Exp.no[1] <- Exp.no.inter[1] * q.freq[s]
      if (j > 1) {
        for (k in 1:(j-1)) {
          Exp.no[k+1] <- Exp.no.inter[k]   * p.freq[s] +
            Exp.no.inter[k+1] * q.freq[s]
        }
      }
      Exp.no[j+1] <- Exp.no.inter[j] * p.freq[s]
      for (k in 1:(j+1)) Exp.no.inter[k] <- Exp.no[k]
    }
    Exp.no
  }

  Compare.obs.exp.freq <- function(Obs.freq, Exp.freq) {
    No <- length(Obs.freq) - 1
    Cum.obs <- array(0, No + 1, dimnames = list(seq(0, No, 1)))
    Cum.exp <- array(0, No + 1, dimnames = list(seq(0, No, 1)))
    Cum.obs[No + 1] <- Obs.freq[No + 1]
    Cum.exp[No + 1] <- Exp.freq[No + 1]
    m <- 1
    m.found <- -1
    if (Cum.obs[No + 1] > Cum.exp[No + 1]) {
      m <- No
      m.found <- 0
    }
    for (j in seq(No, 1, by = -1)) {
      Cum.obs[j] <- Cum.obs[j+1] + Obs.freq[j]
      if (m.found == -1 && Cum.obs[j] > 0) m.found <- 0
      Cum.exp[j] <- Cum.exp[j+1] + Exp.freq[j]
      if (j > 1 && m.found == 0 && Cum.exp[j] > Cum.obs[j]) {
        m <- j
        m.found <- 1
      }
    }
    m
  }

  ## -------- Allocate outputs --------
  Cluster.species <- array(0L, c(n - 1L, n), dimnames = list(as.character(1:(n-1)), species))
  Cluster.info    <- array(0L, c(n - 1L, 2L), dimnames = list(as.character(1:(n-1)), c("k","m")))
  Plot.cluster    <- array(0L, c(N, n - 1L), dimnames = list(as.character(1:N), as.character(1:(n-1))))
  Cluster.merged  <- array(0L, c(n - 1L, 2L))
  Cluster.height  <- numeric(n - 1L)

  vegmatrix2 <- vm
  n2 <- n
  i  <- 0L
  multiple.max <- 1L

  ## -------- Progress bar --------
  pb <- if (progress) utils::txtProgressBar(min = 0, max = n - 1, style = 3) else NULL
  on.exit({ if (!is.null(pb)) close(pb) }, add = TRUE)

  ## -------- Clustering loop --------
  while (i <= (n - 2L)) {
    # full phi similarity among species / clusters
    phi.index <- vegan::designdist(t(vegmatrix2),
                                   method = "(a*d-b*c)/sqrt((a+c)*(b+d)*(a+b)*(c+d))",
                                   terms = c("binary"), alphagamma = FALSE, abcd = TRUE, "phi")
    phi.index[is.na(phi.index)] <- 0
    allmax <- which(phi.index == max(phi.index))

    elements.per.col <- array(0, n2 - 1L)
    for (x in 1:(n2 - 1L)) elements.per.col[x] <- (x - 1L) * n2 - ((x - 1L) * x / 2)
    multiple.max <- length(allmax)
    e1 <- array(0L, multiple.max)
    e2 <- array(0L, multiple.max)
    for (j in 1:multiple.max) {
      e1[j] <- max(which(elements.per.col < allmax[j]))
      e2[j] <- allmax[j] - elements.per.col[e1[j]] + e1[j]
    }

    # de-circularize simultaneous merges
    compare1 <- as.vector(t(as.matrix(cbind(e1, e2))))
    if (anyDuplicated(compare1) > 2) {
      tobekept <- rep(TRUE, multiple.max)
      for (j in 2:multiple.max) {
        compare2 <- compare1[1:(2 * (j - 1))]
        k1 <- length(match(compare2, e1[j])[!is.na(match(compare2, e1[j]))])
        k2 <- length(match(compare2, e2[j])[!is.na(match(compare2, e2[j]))])
        if (k1 > 0 & k2 > 0) tobekept[j] <- FALSE
      }
      e1 <- e1[tobekept]; e2 <- e2[tobekept]; multiple.max <- length(e1)
    }

    i1 <- i + 1L
    for (j in 1:multiple.max) {
      i <- i + 1L
      Cluster.height[i] <- phi.index[allmax][j]

      # left element
      if (substr(colnames(vegmatrix2)[e1[j]], 1, 2) == "c_") {
        cluster.no.1 <- as.numeric(strsplit(colnames(vegmatrix2)[e1[j]], "_")[[1]][2])
        Cluster.merged[i, 1] <- cluster.no.1
        Cluster.species[i, which(Cluster.species[cluster.no.1, ] == 1L)] <- 1L
      } else {
        f1 <- which(colnames(vm) == colnames(vegmatrix2)[e1[j]])
        Cluster.merged[i, 1] <- -f1
        Cluster.species[i, f1] <- 1L
      }

      # right element
      if (substr(colnames(vegmatrix2)[e2[j]], 1, 2) == "c_") {
        cluster.no.2 <- as.numeric(strsplit(colnames(vegmatrix2)[e2[j]], "_")[[1]][2])
        Cluster.merged[i, 2] <- cluster.no.2
        Cluster.species[i, which(Cluster.species[cluster.no.2, ] == 1L)] <- 1L
      } else {
        f2 <- which(colnames(vm) == colnames(vegmatrix2)[e2[j]])
        Cluster.merged[i, 2] <- -f2
        Cluster.species[i, f2] <- 1L
      }

      # rename fused columns to new cluster "c_i"
      cn <- colnames(vegmatrix2)
      cn[c(e1, e2)][cn[c(e1, e2)] == cn[e1[j]]] <- paste0("c_", i)
      cn[c(e1, e2)][cn[c(e1, e2)] == cn[e2[j]]] <- paste0("c_", i)
      colnames(vegmatrix2) <- cn

      # cluster size
      Cluster.info[i, 1] <- sum(Cluster.species[i, ])

      # observed plot frequency
      spp_cols <- which(Cluster.species[i, ] > 0L)
      species.in.plot <- rowSums(vm[, spp_cols, drop = FALSE])
      Obs.plot.freq <- matrix(0, Cluster.info[i, 1] + 1L, dimnames = list(0:Cluster.info[i, 1]))
      tab <- table(species.in.plot)
      Obs.plot.freq[, 1] <- as.numeric(tab[match(rownames(Obs.plot.freq), names(tab))])
      Obs.plot.freq[is.na(Obs.plot.freq)] <- 0

      # expected plot frequency (x N to get counts)
      Exp.plot.freq <- Expected.plot.freq(spp_cols) * N

      # m threshold
      Cluster.info[i, 2] <- Compare.obs.exp.freq(Obs.plot.freq, Exp.plot.freq)

      # plot assignment
      Plot.cluster[species.in.plot >= Cluster.info[i, 2], i] <- 1L
    }

    if (n2 == 3L) {
      name.last.cluster <- colnames(vegmatrix2)[-c(e1[1], e2[1])]
    }

    index.e1 <- !duplicated(colnames(vegmatrix2)[e1], fromLast = TRUE)
    index.e2 <- !duplicated(colnames(vegmatrix2)[e2], fromLast = TRUE)
    index.e  <- c(i1:i)[index.e1 & index.e2]

    vegmatrix2 <- cbind(vegmatrix2[, -unique(c(e1, e2)), drop = FALSE],
                        Plot.cluster[, index.e, drop = FALSE])
    n2 <- ncol(vegmatrix2)
    for (j in seq_along(index.e)) {
      colnames(vegmatrix2)[n2 - length(index.e) + j] <- paste0("c_", index.e[j])
    }
    if (n2 == 2L) colnames(vegmatrix2)[1] <- name.last.cluster

    if (sum(Plot.cluster[, i]) == N) {
      # finalize remaining merges at height 0
      for (j in (i + 1):(n - 1L)) {
        g1 <- ncol(vegmatrix2)
        cluster.no.1 <- as.numeric(strsplit(colnames(vegmatrix2)[g1], "_")[[1]][2])
        Cluster.merged[j, 1] <- cluster.no.1
        g2 <- 1
        cluster.no.2 <- as.numeric(strsplit(colnames(vegmatrix2)[g2], "_")[[1]][2])
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

    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    if (i >= (n - 1L)) break
  }

  ## -------- Result list --------
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

  ## -------- Optional saving --------
  if (!is.null(save_dir) && (save_csv || save_rds)) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (isTRUE(save_csv)) {
      utils::write.csv(Cluster.species, file.path(save_dir, "Cluster.species.csv"), row.names = FALSE)
      utils::write.csv(Cluster.info,    file.path(save_dir, "Cluster.info.csv"),    row.names = FALSE)
      utils::write.csv(Plot.cluster,    file.path(save_dir, "Plot.cluster.csv"),    row.names = FALSE)
      utils::write.csv(Cluster.merged,  file.path(save_dir, "Cluster.merged.csv"),  row.names = FALSE)
    }
    if (isTRUE(save_rds)) {
      # sanitize rds_compress
      valid <- c("xz","gzip","bzip2","none")
      comp <- tolower(rds_compress[1])
      if (!comp %in% valid) {
        warning("rds_compress must be one of ", paste(valid, collapse=", "),
                "; using 'none'.", call. = FALSE)
        comp <- "none"
      }
      saveRDS(res, file = file.path(save_dir, rds_filename), compress = comp)
    }
  }

  res
}
