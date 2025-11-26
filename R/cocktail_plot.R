cocktail_plot <- function(
    x, file = NULL,
    phi_cut        = NULL,
    label_clusters = FALSE,
    cex_species    = 0.3,
    cex_clusters   = 0.6,
    circle_inch    = 0.1,
    page_size      = 300,
    width_in       = NULL,
    height_in      = 10,
    alpha_fill     = 0.18,
    alpha_border   = 0.65,
    palette        = "rainbow",
    png_res        = 300
) {
  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  n  <- ncol(CS)

  species_names <- if (!is.null(x$species)) x$species else colnames(CS)

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
  if (!is.null(phi_cut)) {
    idx <- which(H >= phi_cut)
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
            node_id = i
          )
        }
        bands <- do.call(rbind, blist)
      }
    }
  }

  compute_phi_crossings <- function(ycut, x_from, x_to) {
    if (is.null(bands) || nrow(bands) == 0) return(NULL)
    crosses <- list()

    L_hit <- which(pmin(Pos[, "ly0"], Pos[, "ly1"]) <= ycut &
                     pmax(Pos[, "ly0"], Pos[, "ly1"]) >= ycut)
    if (length(L_hit)) {
      for (r in L_hit) {
        x_leg <- Pos[r, "lx"]
        if (x_leg < x_from || x_leg > x_to) next
        bcand <- which(bands$x0 <= x_leg & bands$x1 >= x_leg)
        if (!length(bcand)) next
        nid <- bands$node_id[bcand[1]]
        crosses[[length(crosses) + 1]] <- data.frame(x = x_leg, y = ycut, label = nid)
      }
    }

    R_hit <- which(pmin(Pos[, "ry0"], Pos[, "ry1"]) <= ycut &
                     pmax(Pos[, "ry0"], Pos[, "ry1"]) >= ycut)
    if (length(R_hit)) {
      for (r in R_hit) {
        x_leg <- Pos[r, "rx"]
        if (x_leg < x_from || x_leg > x_to) next
        bcand <- which(bands$x0 <= x_leg & bands$x1 >= x_leg)
        if (!length(bcand)) next
        nid <- bands$node_id[bcand[1]]
        crosses[[length(crosses) + 1]] <- data.frame(x = x_leg, y = ycut, label = nid)
      }
    }

    if (!length(crosses)) return(NULL)
    do.call(rbind, crosses)
  }

  ## --- pagination ---
  pages <- split(seq_len(n), ceiling(seq_len(n) / page_size))
  if (is.null(width_in)) {
    width_in <- min(150, max(16, 0.06 * max(lengths(pages))))
  }

  draw_page <- function(x_from, x_to) {
    graphics::par(mar = c(8, 5, 1, 1), xaxs = "i", yaxs = "i")
    plot(c(x_from, x_to), c(0.1, -1), type = "n", xaxt = "n", yaxt = "n",
         xlab = "", ylab = expression(paste(phi, " coefficient")))
    graphics::axis(2, las = 2, at = seq(0, -1, by = -0.2), labels = seq(0, 1, by = 0.2))

    if (!is.null(bands) && nrow(bands) > 0) {
      hit <- which(bands$x1 >= x_from & bands$x0 <= x_to)
      for (b in hit) {
        xl <- max(bands$x0[b], x_from) - 0.5
        xr <- min(bands$x1[b], x_to)   + 0.5
        graphics::rect(xl, bands$y0[b], xr, bands$y1[b],
                       col = bands$col[b], border = bands$border[b])
      }
    }

    for (r in 1:nrow(Pos)) {
      lx <- Pos[r, "lx"]; rx <- Pos[r, "rx"]
      if (lx >= x_from && lx <= x_to)
        graphics::segments(lx, Pos[r, "ly0"], lx, Pos[r, "ly1"])
      if (rx >= x_from && rx <= x_to)
        graphics::segments(rx, Pos[r, "ry0"], rx, Pos[r, "ry1"])
    }

    eps <- 1e-9
    hitH <- which(x1 >= x_from - eps & x0 <= x_to + eps)
    for (r in hitH) {
      hx0 <- max(x0[r], x_from)
      hx1 <- min(x1[r], x_to)
      if (hx1 >= hx0 - eps)
        graphics::segments(hx0, y[r], hx1, y[r])
    }

    if (!is.null(phi_cut)) {
      ycuts <- -phi_cut
      ycuts <- ycuts[is.finite(ycuts) & ycuts <= 0 & ycuts >= -1]
      if (length(ycuts)) {
        for (yy in ycuts) {
          graphics::segments(x_from - 0.5, yy, x_to + 0.5, yy,
                             lty = 2, lwd = 1, col = "grey20")
        }
      }
    }

    idx <- x_from:x_to
    graphics::axis(1, at = idx, labels = species_names[species_order][idx],
                   las = 2, cex.axis = cex_species)

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
    } else if (identical(label_clusters, "phi_cut")) {
      if (!is.null(phi_cut)) {
        ycut <- -phi_cut
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

  starts  <- seq(1, n, by = page_size)
  ends    <- pmin(starts + page_size - 1, n)
  n_pages <- length(starts)

  use_current_dev <- FALSE
  if (is.null(file)) {
    use_current_dev <- TRUE
  } else {
    is_pdf <- grepl("\\.pdf$", file, ignore.case = TRUE)
    is_png <- grepl("\\.png$", file, ignore.case = TRUE)
    if (!is_pdf && !is_png) {
      warning(
        "`file` does not end with '.pdf' or '.png'. ",
        "Plot was drawn on the current device instead. ",
        "To save the plot to disk, supply a filename ending in '.pdf' or '.png'."
      )
      use_current_dev <- TRUE
    }
  }

  if (use_current_dev) {
    p <- 1L
    draw_page(starts[p], ends[p])
    graphics::mtext(sprintf("Species %d-%d of %d", starts[p], ends[p], n),
                    side = 3, line = 0.2, cex = 0.8)
    return(invisible(NULL))
  }

  if (grepl("\\.pdf$", file, ignore.case = TRUE)) {
    grDevices::pdf(file, width = width_in, height = height_in, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)

    for (p in seq_along(starts)) {
      draw_page(starts[p], ends[p])
      graphics::mtext(sprintf("Species %d-%d of %d", starts[p], ends[p], n),
                      side = 3, line = 0.2, cex = 0.8)
    }

  } else if (grepl("\\.png$", file, ignore.case = TRUE)) {
    base <- sub("\\.png$", "", file, ignore.case = TRUE)

    for (p in seq_along(starts)) {
      fname <- if (n_pages == 1) {
        file
      } else {
        sprintf("%s_page%02d.png", base, p)
      }

      grDevices::png(
        filename = fname,
        width  = width_in,
        height = height_in,
        units  = "in",
        res    = png_res
      )

      draw_page(starts[p], ends[p])
      graphics::mtext(sprintf("Species %d-%d of %d", starts[p], ends[p], n),
                      side = 3, line = 0.2, cex = 0.8)

      grDevices::dev.off()
    }
  }

  invisible(NULL)
}
