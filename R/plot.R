#' Plot Cocktail dendrogram to PDF
#'
#' @param x cocktail object (from \code{cocktail_cluster()})
#' @param file output PDF path
#' @param page_size tips per page (use \code{ncol(x$Cluster.species)} for one page)
#' @param width_in PDF width in inches (NULL = auto; capped at 150")
#' @param height_in PDF height in inches
#' @param cex_labels axis label size
#' @param bands_phi numeric or NULL. If given, draw background bands for
#'  "parent" clusters with \eqn{\phi \ge bands_phi}.
#' @param palette palette name for \code{hcl.colors()} (e.g. "Set3","Dark2","Paired") or "rainbow"
#' @param alpha_fill,alpha_border transparency for fill/border
#'
#' @importFrom graphics axis par segments title rect mtext
#' @importFrom stats setNames
#' @export
plot_cocktail <- function(
    x, file,
    page_size  = 300,
    width_in   = NULL,
    height_in  = 10,
    cex_labels = 0.3,
    bands_phi  = NULL,
    palette    = "Set3",
    alpha_fill   = 0.18,
    alpha_border = 0.65
) {
  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  n  <- ncol(CS)
  species_names <- colnames(CS)

  ## --- species order (as in original code) ---
  ord <- order(H)
  Species.sort <- vapply(seq_len(n), function(i) paste(CS[ord, i], collapse = ""), "")
  species_order <- order(Species.sort)

  ## --- precompute coordinates for every merge (legs + y) ---
  Pos <- matrix(NA_real_, nrow = n - 1, ncol = 6,
                dimnames = list(1:(n - 1), c("lx","ly0","ly1","rx","ry0","ry1")))
  for (i in 1:(n - 1)) {
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
  y  <- -H[seq_len(n - 1)]

  ## --- build bands from phi cut (parent clusters only) ---
  bands <- NULL
  if (!is.null(bands_phi)) {
    idx <- which(H >= bands_phi)
    if (length(idx)) {
      # Drop any node that is a child of another above-threshold node
      children <- unique(as.integer(CM[idx, ][CM[idx, ] > 0]))
      top_idx  <- sort(setdiff(idx, intersect(idx, children)))
      if (length(top_idx)) {
        k <- length(top_idx)
        cols_base <-
          if (palette %in% rownames(grDevices::hcl.pals())) {
            grDevices::hcl.colors(max(k, 3), palette = palette)
          } else if (tolower(palette) == "rainbow") {
            grDevices::rainbow(k)
          } else {
            grDevices::hcl.colors(max(k, 3), palette = "Set3")
          }
        cols_base <- cols_base[seq_len(k)]

        blist <- vector("list", k)
        for (ii in seq_along(top_idx)) {
          i <- top_idx[ii]
          spp  <- which(CS[i, ] == 1)
          tips <- sort(match(spp, species_order))
          if (length(tips) < 2) next
          blist[[ii]] <- data.frame(
            x0 = min(tips), x1 = max(tips),
            y0 = -1, y1 = 0.1,
            col    = grDevices::adjustcolor(cols_base[ii], alpha.f = alpha_fill),
            border = grDevices::adjustcolor(cols_base[ii], alpha.f = alpha_border)
          )
        }
        bands <- do.call(rbind, blist)
      }
    }
  }

  ## --- pagination + device ---
  pages <- split(seq_len(n), ceiling(seq_len(n) / page_size))
  if (is.null(width_in)) {
    width_in <- min(150, max(16, 0.06 * max(lengths(pages))))  # cap wide PDFs
  }
  grDevices::pdf(file, width = width_in, height = height_in, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  draw_page <- function(x_from, x_to) {
    graphics::par(mar = c(8, 5, 1, 1), xaxs = "i", yaxs = "i")
    plot(c(x_from, x_to), c(0.1, -1), type = "n", xaxt = "n", yaxt = "n",
         xlab = "", ylab = expression(paste(phi, " coefficient")))
    graphics::axis(2, las = 2, at = seq(0, -1, by = -0.2), labels = seq(0, 1, by = 0.2))

    # background bands (clip to page; +/- 0.5 to cover tick centers)
    if (!is.null(bands) && nrow(bands) > 0) {
      hit <- which(bands$x1 >= x_from & bands$x0 <= x_to)
      for (b in hit) {
        xl <- max(bands$x0[b], x_from) - 0.5
        xr <- min(bands$x1[b], x_to)   + 0.5
        graphics::rect(xl, bands$y0[b], xr, bands$y1[b],
                       col = bands$col[b], border = bands$border[b])
      }
    }

    # vertical legs (only those inside window)
    for (r in 1:nrow(Pos)) {
      lx <- Pos[r, "lx"]; rx <- Pos[r, "rx"]
      if (lx >= x_from && lx <= x_to) graphics::segments(lx, Pos[r, "ly0"], lx, Pos[r, "ly1"])
      if (rx >= x_from && rx <= x_to) graphics::segments(rx, Pos[r, "ry0"], rx, Pos[r, "ry1"])
    }

    # horizontal connectors with robust clipping
    eps <- 1e-9
    hitH <- which(x1 >= x_from - eps & x0 <= x_to + eps)
    for (r in hitH) {
      hx0 <- max(x0[r], x_from)
      hx1 <- min(x1[r], x_to)
      if (hx1 >= hx0 - eps) graphics::segments(hx0, y[r], hx1, y[r])
    }

    # tip labels
    idx <- x_from:x_to
    graphics::axis(1, at = idx, labels = species_names[species_order][idx],
                   las = 2, cex.axis = cex_labels)
  }

  starts <- seq(1, n, by = page_size)
  ends   <- pmin(starts + page_size - 1, n)
  for (p in seq_along(starts)) {
    draw_page(starts[p], ends[p])
    graphics::mtext(sprintf("Species %d-%d of %d", starts[p], ends[p], n),
                    side = 3, line = 0.2, cex = 0.8)
  }

  # optional: write species-order CSV next to the PDF
  ord_tbl <- data.frame(x = seq_len(n), species = species_names[species_order])
  utils::write.csv(ord_tbl, sub("\\.pdf$", "_species_order.csv", file), row.names = FALSE)
}
