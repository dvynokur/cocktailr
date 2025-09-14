#' Fuzzy Cocktail memberships (no thresholding)
#'
#' @param x Result of `cocktail_cluster()`.
#' @param vegmatrix plots x species cover (numeric); used to get p/a for phi.
#' @param cuts numeric vector of phi cuts that define the parent groups (e.g., c(0.2,0.25,0.3)).
#' @return data.frame with columns: cut, group, species, phi
#' @export
cocktail_fuzzy <- function(x, vegmatrix, cuts = 0.30) {
  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  species <- colnames(CS)

  # species order (same as plotting) - needed to define groups consistently
  ord <- order(H)
  Species.sort <- vapply(seq_len(ncol(CS)), function(i) paste(CS[ord, i], collapse=""), "")
  species_order <- order(Species.sort)

  # p/a matrix for phi calculation (plots x species)
  vm <- as.matrix(vegmatrix)
  if (is.null(colnames(vm))) stop("vegmatrix must have species column names.")
  vm_bin <- (vm > 0) * 1L
  # keep only species in CS
  common <- intersect(colnames(vm_bin), species)
  if (!length(common)) stop("No overlapping species between vegmatrix and clustering result.")
  vm_bin <- vm_bin[, common, drop = FALSE]

  # fast phi for one species vs a set of plots
  .phi_one <- function(pa, in_grp) {
    # pa: vector (plots) 0/1 for a species
    # in_grp: logical vector (plots) TRUE if plot %in% group
    a <- sum(pa == 1 &  in_grp)
    b <- sum(pa == 1 & !in_grp)
    c <- sum(pa == 0 &  in_grp)
    d <- sum(pa == 0 & !in_grp)
    den <- sqrt((a + c) * (b + d) * (a + b) * (c + d))
    if (den == 0) return(0)
    (a * d - b * c) / den
  }

  out_list <- list()
  gi <- 0

  for (cut in cuts) {
    # find "topmost" clusters with height >= cut
    idx <- which(H >= cut)
    if (!length(idx)) next
    children <- unique(as.integer(CM[idx, ][CM[idx, ] > 0]))
    top_idx <- sort(setdiff(idx, intersect(idx, children)))
    if (!length(top_idx)) next

    # build plot-membership for each cluster: Plot.cluster[,i] == 1
    # (if available in x), else approximate: plots where >= m species present
    if (!is.null(x$Plot.cluster)) {
      PC <- x$Plot.cluster  # N x (n-1)
    } else {
      # slow fallback: compute from CS & vegmatrix if Plot.cluster missing
      stop("x$Plot.cluster not found; cocktail_cluster() should return it.")
    }

    for (ii in seq_along(top_idx)) {
      i <- top_idx[ii]
      grp_plots <- as.logical(PC[, i])
      if (!any(grp_plots)) next

      # compute phi for each species (only those in vm_bin)
      phi_vec <- vapply(colnames(vm_bin), function(sp) .phi_one(vm_bin[, sp], grp_plots), numeric(1))
      gi <- gi + 1L
      out_list[[gi]] <- data.frame(
        cut = cut,
        group = paste0("G", ii, "_cut", formatC(cut, format = "f", digits = 2)),
        species = names(phi_vec),
        phi = as.numeric(phi_vec),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(out_list)) {
    return(data.frame(cut = numeric(0), group = character(0), species = character(0), phi = numeric(0)))
  }
  do.call(rbind, out_list)
}
