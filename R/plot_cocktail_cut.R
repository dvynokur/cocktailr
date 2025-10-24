#' Plot Cocktail dendrogram to PDF
#'
#' @description
#' Draw a dendrogram-style plot of a fitted Cocktail clustering and write it to
#' a PDF file. Species (tips) are arranged left-to-right in the same order used
#' by the original algorithm, and the y-axis shows the \eqn{\phi} height of each
#' merge. Optionally, highlight “parent” clusters at or above a chosen cut and
#' draw a dashed horizontal line at that cut.
#'
#' @param x A Cocktail object (e.g., from \code{cocktail_cluster()}), with
#'   components \code{Cluster.species}, \code{Cluster.merged}, and \code{Cluster.height}.
#' @param file Path to the output PDF.
#' @param page_size Integer, number of tips per page. Use
#'   \code{ncol(x$Cluster.species)} to place all tips on a single page.
#' @param width_in PDF width in inches. If \code{NULL}, a width is chosen
#'   automatically (capped at 150 inches) based on \code{page_size}.
#' @param height_in PDF height in inches (default 10).
#' @param cex_labels Numeric expansion factor for tip labels on the x-axis.
#' @param bands_phi Numeric or \code{NULL}. If given, draw background bands for
#'   parent clusters with \eqn{\phi \ge \mathrm{bands\_phi}}, and also draw a
#'   dashed horizontal line at this \eqn{\phi} level across each page.
#' @param palette Color palette name for \code{grDevices::hcl.colors()}
#'   (e.g., \code{"Set3"}, \code{"Dark2"}, \code{"Paired"}) or \code{"rainbow"}.
#' @param alpha_fill,alpha_border Numeric alpha levels (0–1) for band fill and border.
#' @param label_clusters Logical/character. One of \code{FALSE} (default),
#'   \code{"all"}, or \code{"phi"}. If \code{"all"}, label every internal node
#'   (merge). If \code{"phi"}, label cluster numbers at the intersections of the
#'   cut \code{bands_phi} with vertical branches (parent clusters at/above
#'   \code{bands_phi}). Labels correspond to node indices (same as \code{assign_releves()}),
#'   but shown without the \code{"c_"} prefix.
#' @param cex_clusters Numeric expansion factor for cluster-number text (default 0.6).
#' @param circle_inch Circle radius (inches) for the white/black bubble (default 0.1).
#'
#' @details
#' Species order is reproduced from the original Cocktail implementation by
#' sorting the per-tip membership strings derived from \code{x$Cluster.species}.
#' Large trees are automatically paginated using \code{page_size}; the device
#' width grows with the number of tips on a page but is capped to avoid
#' unwieldy PDFs.
#'
#' @return Creates a PDF on disk; returns \code{NULL} invisibly.
#'
#' @importFrom graphics axis par segments title rect mtext symbols text
#' @importFrom stats setNames
#' @export
plot_cocktail_cut <- function(
    x, file,
    page_size    = 300,
    width_in     = NULL,
    height_in    = 10,
    cex_labels   = 0.3,
    bands_phi    = NULL,
    palette      = "Set3",
    alpha_fill   = 0.18,
    alpha_border = 0.65,
    label_clusters = FALSE,
    cex_clusters  = 0.6,
    circle_inch   = 0.1
) {
  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  n  <- ncol(CS)
  species_names <- colnames(CS)

  ## --- species order ---
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
  top_idx <- integer(0)
  if (!is.null(bands_phi)) {
    idx <- which(H >= bands_phi)
    if (length(idx)) {
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
            border = grDevices::adjustcolor(cols_base[ii], alpha.f = alpha_border),
            node_id = i       # use actual node index
          )
        }
        bands <- do.call(rbind, blist)
      }
    }
  }

  ## --- helper: compute cluster labels at phi cut ---
  compute_phi_crossings <- function(ycut, x_from, x_to) {
    if (is.null(bands) || nrow(bands) == 0) return(NULL)
    crosses <- list()

    L_hit <- which(pmin(Pos[, "ly0"], Pos[, "ly1"]) <= ycut & pmax(Pos[, "ly0"], Pos[, "ly1"]) >= ycut)
    if (length(L_hit)) {
      for (r in L_hit) {
        x_leg <- Pos[r, "lx"]
        if (x_leg < x_from || x_leg > x_to) next
        bcand <- which(bands$x0 <= x_leg & bands$x1 >= x_leg)
        if (!length(bcand)) next
        nid <- bands$node_id[bcand[1]]   # node index
        crosses[[length(crosses) + 1]] <- data.frame(x = x_leg, y = ycut, label = nid)
      }
    }

    R_hit <- which(pmin(Pos[, "ry0"], Pos[, "ry1"]) <= ycut & pmax(Pos[, "ry0"], Pos[, "ry1"]) >= ycut)
    if (length(R_hit)) {
      for (r in R_hit) {
        x_leg <- Pos[r, "rx"]
        if (x_leg < x_from || x_leg > x_to) next
        bcand <- which(bands$x0 <= x_leg & bands$x1 >= x_leg)
        if (!length(bcand)) next
        nid <- bands$node_id[bcand[1]]   # node index
        crosses[[length(crosses) + 1]] <- data.frame(x = x_leg, y = ycut, label = nid)
      }
    }

    if (!length(crosses)) return(NULL)
    do.call(rbind, crosses)
  }

  ## --- pagination + device ---
  pages <- split(seq_len(n), ceiling(seq_len(n) / page_size))
  if (is.null(width_in)) {
    width_in <- min(150, max(16, 0.06 * max(lengths(pages))))
  }
  grDevices::pdf(file, width = width_in, height = height_in, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  draw_page <- function(x_from, x_to) {
    graphics::par(mar = c(8, 5, 1, 1), xaxs = "i", yaxs = "i")
    plot(c(x_from, x_to), c(0.1, -1), type = "n", xaxt = "n", yaxt = "n",
         xlab = "", ylab = expression(paste(phi, " coefficient")))
    graphics::axis(2, las = 2, at = seq(0, -1, by = -0.2), labels = seq(0, 1, by = 0.2))

    # background bands
    if (!is.null(bands) && nrow(bands) > 0) {
      hit <- which(bands$x1 >= x_from & bands$x0 <= x_to)
      for (b in hit) {
        xl <- max(bands$x0[b], x_from) - 0.5
        xr <- min(bands$x1[b], x_to)   + 0.5
        graphics::rect(xl, bands$y0[b], xr, bands$y1[b],
                       col = bands$col[b], border = bands$border[b])
      }
    }

    # vertical legs
    for (r in 1:nrow(Pos)) {
      lx <- Pos[r, "lx"]; rx <- Pos[r, "rx"]
      if (lx >= x_from && lx <= x_to) graphics::segments(lx, Pos[r, "ly0"], lx, Pos[r, "ly1"])
      if (rx >= x_from && rx <= x_to) graphics::segments(rx, Pos[r, "ry0"], rx, Pos[r, "ry1"])
    }

    # horizontal connectors
    eps <- 1e-9
    hitH <- which(x1 >= x_from - eps & x0 <= x_to + eps)
    for (r in hitH) {
      hx0 <- max(x0[r], x_from)
      hx1 <- min(x1[r], x_to)
      if (hx1 >= hx0 - eps) graphics::segments(hx0, y[r], hx1, y[r])
    }

    # dashed cut line(s)
    if (!is.null(bands_phi)) {
      ycuts <- -bands_phi
      ycuts <- ycuts[is.finite(ycuts) & ycuts <= 0 & ycuts >= -1]
      if (length(ycuts)) {
        for (yy in ycuts) {
          graphics::segments(x_from - 0.5, yy, x_to + 0.5, yy, lty = 2, lwd = 1, col = "grey20")
        }
      }
    }

    # tip labels
    idx <- x_from:x_to
    graphics::axis(1, at = idx, labels = species_names[species_order][idx],
                   las = 2, cex.axis = cex_labels)

    ## --- CLUSTER / NODE LABELS ---
    if (identical(label_clusters, "all")) {
      hit_nodes <- which(x1 >= x_from & x0 <= x_to)
      if (length(hit_nodes)) {
        xm <- (pmax(x0[hit_nodes], x_from) + pmin(x1[hit_nodes], x_to)) / 2
        ym <- y[hit_nodes]
        graphics::symbols(xm, ym,
                          circles = rep(1, length(xm)),
                          inches  = circle_inch,
                          add = TRUE, bg = "white", fg = "black", lwd = 0.6)
        graphics::text(xm, ym, labels = hit_nodes, cex = cex_clusters)
      }
    } else if (identical(label_clusters, "phi")) {
      if (!is.null(bands_phi) && length(top_idx)) {
        ycut <- -bands_phi
        for (yy in ycut) {
          cross <- compute_phi_crossings(yy, x_from, x_to)
          if (!is.null(cross) && nrow(cross)) {
            graphics::symbols(cross$x, cross$y,
                              circles = rep(1, nrow(cross)),
                              inches  = circle_inch,
                              add = TRUE, bg = "white", fg = "black", lwd = 0.6)
            graphics::text(cross$x, cross$y,
                           labels = cross$label,
                           cex = cex_clusters)
          }
        }
      }
    }
  }

  starts <- seq(1, n, by = page_size)
  ends   <- pmin(starts + page_size - 1, n)
  for (p in seq_along(starts)) {
    draw_page(starts[p], ends[p])
    graphics::mtext(sprintf("Species %d-%d of %d", starts[p], ends[p], n),
                    side = 3, line = 0.2, cex = 0.8)
  }
}
