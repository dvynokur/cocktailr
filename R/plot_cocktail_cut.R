#' Plot Cocktail dendrogram from a phi cut (experimental)
#'
#' @inheritParams plot_cocktail
#' @param phi_cut Numeric in (-1,1); draw only merges at/below this height and
#'   label tips by the parent groups at/above this cut.
#' @param tip_label "c_id" (default) or "id" for group labels at tips.
#' @param show_k Logical; if TRUE, append " (k=<n>)" to parent labels.
#' @param show_singletons Logical; if TRUE (default) show species that are not in
#'   any parent at the cut as grey tiles (without trunks) and label them.
#'
#' @return Invisibly returns NULL after writing a PDF.
#' @export
plot_cocktail_cut <- function(
    x, file,
    phi_cut,
    tip_label       = c("c_id","id"),
    show_k          = TRUE,
    show_singletons = TRUE,
    width_in   = NULL,
    height_in  = 10,
    cex_labels = 0.7,
    palette    = "Set3",
    alpha_fill   = 0.18,
    alpha_border = 0.65
) {
  tip_label <- match.arg(tip_label)

  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  n  <- ncol(CS)
  sp_names <- colnames(CS); if (is.null(sp_names)) sp_names <- paste0("sp_", seq_len(n))

  ## allow negatives now
  if (!is.numeric(phi_cut) || length(phi_cut) != 1L || !is.finite(phi_cut) ||
      phi_cut < -1 || phi_cut > 1) stop("`phi_cut` must be in [-1, 1].")

  ## --- species order (as in original) ---
  ord <- order(H)
  Species.sort  <- vapply(seq_len(n), function(i) paste(CS[ord, i], collapse=""), "")
  species_order <- order(Species.sort)

  ## --- original coordinates for every merge (same as plot_cocktail) ---
  m <- nrow(CM)
  Pos <- matrix(NA_real_, nrow = m, ncol = 6,
                dimnames = list(1:m, c("lx","ly0","ly1","rx","ry0","ry1")))
  for (i in 1:m) {
    if (CM[i, 1] < 0) {
      Pos[i, "lx"]  <- which(species_order == -CM[i, 1])
      Pos[i, "ly0"] <- -1
      Pos[i, "ly1"] <- -H[i]
    } else {
      spp <- which(CS[CM[i, 1], ] == 1)
      Pos[i, "lx"]  <- mean(match(spp, species_order))
      Pos[i, "ly0"] <- -H[CM[i, 1]]
      Pos[i, "ly1"] <- -H[i]
    }
    if (CM[i, 2] < 0) {
      Pos[i, "rx"]  <- which(species_order == -CM[i, 2])
      Pos[i, "ry0"] <- -1
      Pos[i, "ry1"] <- -H[i]
    } else {
      spp <- which(CS[CM[i, 2], ] == 1)
      Pos[i, "rx"]  <- mean(match(spp, species_order))
      Pos[i, "ry0"] <- -H[CM[i, 2]]
      Pos[i, "ry1"] <- -H[i]
    }
  }
  x0 <- pmin(Pos[, "lx"], Pos[, "rx"])
  x1 <- pmax(Pos[, "lx"], Pos[, "rx"])
  y  <- -H[seq_len(m)]

  ## --- parent groups at/above the cut (top-most only) ---
  idx <- which(H >= phi_cut)
  kids <- if (length(idx)) unique(as.integer(CM[idx, ][CM[idx, ] > 0])) else integer(0)
  parents <- sort(setdiff(idx, intersect(idx, kids)))
  parent_labels <- paste0("c_", parents)

  # map species -> parent id (NA if none)
  parent_map <- rep(NA_integer_, n)
  for (p in parents) parent_map[CS[p, ] == 1L] <- p

  # parent spans
  parent_spans <- vapply(parents, function(p) {
    tips <- sort(match(which(CS[p, ] == 1L), species_order))
    c(min(tips), max(tips))
  }, numeric(2))
  parent_spans <- t(parent_spans)
  colnames(parent_spans) <- c("min","max")
  rownames(parent_spans) <- parents

  # labels for parents
  labs <- if (tip_label == "c_id") paste0("c_", parents) else as.character(parents)
  if (isTRUE(show_k)) {
    kvec <- rowSums(CS[parents, , drop = FALSE] == 1L)
    labs <- paste0(labs, " (k=", kvec, ")")
  }

  # optional singletons
  single_idx <- which(!is.finite(parent_map))
  single_centers <- if (length(single_idx)) match(single_idx, species_order) else integer(0)

  ## --- device / window ---
  n_tips <- length(parents) + if (show_singletons) length(single_idx) else 0L
  if (is.null(width_in)) width_in <- min(150, max(10, 0.12 * max(n_tips, 1)))
  grDevices::pdf(file, width = width_in, height = height_in, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(mar = c(6, 5, 2, 2), xaxs = "i", yaxs = "i")

  ## allow negative depths: bottom extends to the minimum φ seen (or 0 if all >0)
  min_phi <- min(H, 0, na.rm = TRUE)
  xlim <- c(0.5, n + 0.5)
  ylim <- c(-phi_cut, -min_phi)

  plot(xlim, ylim, type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = expression(paste(phi, " coefficient")))
  yt <- pretty(c(min_phi, phi_cut))
  graphics::axis(2, las = 2, at = -yt, labels = yt)

  ## --- colored tiles for parent spans + grey tiles for singletons -----------
  if (alpha_fill > 0 || alpha_border > 0) {
    k <- max(1L, length(parents))
    cols_base <-
      if (palette %in% rownames(grDevices::hcl.pals())) {
        grDevices::hcl.colors(max(k, 3), palette = palette)
      } else if (tolower(palette) == "rainbow") {
        grDevices::rainbow(k)
      } else {
        grDevices::hcl.colors(max(k, 3), palette = "Set3")
      }
    cols_base <- cols_base[seq_len(k)]
    for (ii in seq_along(parents)) {
      spn <- parent_spans[ii, ]
      graphics::rect(spn["min"] - 0.5, -phi_cut, spn["max"] + 0.5, -phi_cut + 0.08,
                     col    = grDevices::adjustcolor(cols_base[ii], alpha.f = alpha_fill),
                     border = grDevices::adjustcolor(cols_base[ii], alpha.f = alpha_border))
    }
    if (isTRUE(show_singletons) && length(single_centers)) {
      for (cx in single_centers) {
        graphics::rect(cx - 0.5, -phi_cut, cx + 0.5, -phi_cut + 0.08,
                       col    = grDevices::adjustcolor("grey70", alpha.f = alpha_fill),
                       border = grDevices::adjustcolor("grey40", alpha.f = alpha_border))
      }
    }
  }

  ## --- draw original segments, clipped at the cut ---------------------------
  draw_idx <- which(H <= phi_cut)  # merges to show (at/under the cut; includes negatives if asked)
  if (length(draw_idx)) {
    for (r in draw_idx) {
      # vertical legs clipped to the cut
      ly0 <- max(-phi_cut, Pos[r, "ly0"])
      ry0 <- max(-phi_cut, Pos[r, "ry0"])
      graphics::segments(Pos[r, "lx"], ly0, Pos[r, "lx"], Pos[r, "ly1"])
      graphics::segments(Pos[r, "rx"], ry0, Pos[r, "rx"], Pos[r, "ry1"])
      # horizontal connector
      graphics::segments(x0[r], y[r], x1[r], y[r])
    }
  }

  ## --- tip labels (parents + optional singletons) ---------------------------
  centers <- if (length(parents)) vapply(seq_len(nrow(parent_spans)), function(i) mean(parent_spans[i, ]), 0.0) else numeric(0)
  at <- centers; lab <- labs
  if (isTRUE(show_singletons) && length(single_centers)) {
    at  <- c(at, single_centers)
    lab <- c(lab, sp_names[species_order][single_centers])
  }
  if (length(at)) {
    ord_lab <- order(at); at <- at[ord_lab]; lab <- lab[ord_lab]
    graphics::axis(1, at = at, labels = lab, las = 2, cex.axis = cex_labels)
  }

  # dashed cut
  graphics::segments(xlim[1], -phi_cut, xlim[2], -phi_cut, lty = 2, col = "grey30")

  invisible(NULL)
}
