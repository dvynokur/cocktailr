#' Plot Cocktail dendrogram (PDF/PNG or current device)
#'
#' @description
#' Draw a dendrogram-style plot of a fitted Cocktail clustering.
#' Species (tips) are arranged left-to-right in the same order used
#' by the original algorithm, and the y-axis shows the \eqn{\phi}
#' height of each merge.
#'
#' If \code{file} ends with \code{".pdf"} or \code{".png"}, the plot
#' is written to disk. Otherwise (including \code{file = NULL}),
#' the first page is drawn on the current graphics device (e.g.,
#' the Plots panel in RStudio).
#'
#' @param x A Cocktail object (e.g. from \code{cocktail_cluster()}),
#'   with components \code{Cluster.species}, \code{Cluster.merged},
#'   and \code{Cluster.height}. If present, \code{x$species} is used
#'   as the species names; otherwise \code{colnames(Cluster.species)}.
#' @param file Optional path to output. If it ends in:
#'   \itemize{
#'     \item \code{".pdf"} — create a (multi-page) PDF;
#'     \item \code{".png"} — create one or more PNG files
#'           (\code{*_page01.png}, \code{*_page02.png}, … if multiple pages);
#'     \item anything else or \code{NULL} — draw only the first page
#'           on the current graphics device.
#'   }
#' @param phi_cut Optional numeric. If given, draw background bands for
#'   parent clusters with \eqn{\phi \ge \mathrm{phi\_cut}}, and a dashed
#'   horizontal line at this \eqn{\phi} level across each page.
#' @param label_clusters Logical/character. One of:
#'   \itemize{
#'     \item \code{FALSE} — no internal node labels;
#'     \item \code{"all"} — label every internal node, placing the label
#'           where the parent branch meets the node (intersection of the
#'           parent’s vertical leg with the parent’s horizontal level;
#'           for root nodes, at their own merge height), \emph{plus} an
#'           additional set of labels at \eqn{\phi = 0} corresponding to
#'           the \eqn{\phi_\mathrm{cut} = 0} clusters (i.e. the same labels
#'           as \code{label_clusters = "phi_cut"} with \code{phi_cut = 0});
#'     \item \code{"phi_cut"} — only label nodes whose vertical branch
#'           crosses \code{phi_cut} (as in the dashed line).
#'   }
#' @param cex_species Expansion factor for species (tip) labels on the x-axis.
#' @param cex_clusters Expansion factor for cluster-number labels.
#' @param circle_inch Circle radius (inches) for the white/black bubble
#'   behind cluster labels.
#' @param page_size Integer, number of tips per page. Use
#'   \code{ncol(x$Cluster.species)} to place all tips on a single page.
#' @param width_in Plot width in inches. If \code{NULL}, a width is chosen
#'   automatically (capped at 150 inches) based on \code{page_size}.
#' @param height_in Plot height in inches (default 10).
#' @param alpha_fill,alpha_border Numeric alpha levels (0–1) for band fill
#'   and border.
#' @param palette Color palette name for \code{grDevices::hcl.colors()}
#'   (e.g., \code{"Set3"}, \code{"Dark2"}, \code{"Paired"}) or \code{"rainbow"}.
#'   Default \code{"rainbow"}.
#' @param png_res Resolution in dpi for PNG output (default 300).
#'
#' @return Invisibly returns \code{NULL}. Produces a plot on the current
#'   device or writes PDF/PNG files, depending on \code{file}.
#'
#' @importFrom graphics axis par segments mtext rect symbols text
#' @importFrom grDevices hcl.colors adjustcolor pdf png dev.off hcl.pals rainbow
#' @export
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

  ## --- species order (Cocktail-like) ---
  ord <- order(H)
  Species.sort <- vapply(seq_len(n), function(i) paste(CS[ord, i], collapse = ""), "")
  species_order <- order(Species.sort)

  ## --- precompute coordinates for every merge (legs + y) ---
  Pos <- matrix(
    NA_real_,
    nrow = n - 1,
    ncol = 6,
    dimnames = list(1:(n - 1), c("lx","ly0","ly1","rx","ry0","ry1"))
  )

  for (i in 1:(n - 1)) {
    # left child
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

    # right child
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

  ## --- cluster "centre" (horizontal midpoint of tips) ---
  center_x <- (x0 + x1) / 2

  ## --- parent index for each node (1..n-1), NA for roots ---
  parent_idx <- rep(NA_integer_, nrow(CM))
  for (p in seq_len(nrow(CM))) {
    for (j in 1:2) {
      ch <- CM[p, j]
      if (ch > 0L && ch <= nrow(CM)) {
        parent_idx[ch] <- p
      }
    }
  }

  ## --- helper: remove overlapping labels, keep higher ones ---
  filter_non_overlapping <- function(df, rx, ry) {
    if (is.null(df) || nrow(df) <= 1L) return(df)
    # higher (coarser) = y closer to 0 (less negative)
    ord <- order(df$y, decreasing = TRUE)
    df  <- df[ord, , drop = FALSE]

    keep <- logical(nrow(df))
    keep[1L] <- TRUE

    for (i in 2:nrow(df)) {
      xi <- df$x[i]
      yi <- df$y[i]
      overlap <- FALSE
      for (j in which(keep)) {
        if (abs(xi - df$x[j]) < 2 * rx &&
            abs(yi - df$y[j]) < 2 * ry) {
          overlap <- TRUE
          break
        }
      }
      if (!overlap) keep[i] <- TRUE
    }
    df[keep, , drop = FALSE]
  }

  ## --- build bands from phi_cut (parent clusters only) ---
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

  ## --- EXTRA bands for phi = 0 labelling in "all" mode (as if phi_cut = 0) ---
  bands_phi0 <- NULL
  if (identical(label_clusters, "all")) {
    idx0 <- which(H >= 0)
    if (length(idx0)) {
      children0 <- unique(as.integer(CM[idx0, ][CM[idx0, ] > 0]))
      top_idx0  <- sort(setdiff(idx0, intersect(idx0, children0)))

      if (length(top_idx0)) {
        blist0 <- vector("list", length(top_idx0))
        for (ii in seq_along(top_idx0)) {
          i <- top_idx0[ii]
          spp  <- which(CS[i, ] == 1)
          tips <- sort(match(spp, species_order))
          if (length(tips) < 2) next
          blist0[[ii]] <- data.frame(
            x0 = min(tips), x1 = max(tips),
            y0 = -1, y1 = 0.1,
            node_id = i
          )
        }
        bands_phi0 <- do.call(rbind, blist0)
      }
    }
  }

  ## --- helper: cluster labels at phi_cut crossings (vertical branches) ---
  compute_phi_crossings <- function(ycut, x_from, x_to, bands_df) {
    if (is.null(bands_df) || nrow(bands_df) == 0) return(NULL)
    crosses <- list()

    # left legs
    L_hit <- which(
      pmin(Pos[, "ly0"], Pos[, "ly1"]) <= ycut &
        pmax(Pos[, "ly0"], Pos[, "ly1"]) >= ycut
    )
    if (length(L_hit)) {
      for (r in L_hit) {
        x_leg <- Pos[r, "lx"]
        if (x_leg < x_from || x_leg > x_to) next
        bcand <- which(bands_df$x0 <= x_leg & bands_df$x1 >= x_leg)
        if (!length(bcand)) next
        nid <- bands_df$node_id[bcand[1]]
        crosses[[length(crosses) + 1]] <- data.frame(
          x = x_leg, y = ycut, label = nid
        )
      }
    }

    # right legs
    R_hit <- which(
      pmin(Pos[, "ry0"], Pos[, "ry1"]) <= ycut &
        pmax(Pos[, "ry0"], Pos[, "ry1"]) >= ycut
    )
    if (length(R_hit)) {
      for (r in R_hit) {
        x_leg <- Pos[r, "rx"]
        if (x_leg < x_from || x_leg > x_to) next
        bcand <- which(bands_df$x0 <= x_leg & bands_df$x1 >= x_leg)
        if (!length(bcand)) next
        nid <- bands_df$node_id[bcand[1]]
        crosses[[length(crosses) + 1]] <- data.frame(
          x = x_leg, y = ycut, label = nid
        )
      }
    }

    if (!length(crosses)) return(NULL)
    do.call(rbind, crosses)
  }

  ## --- pagination (+ auto width if needed) ---
  pages <- split(seq_len(n), ceiling(seq_len(n) / page_size))
  if (is.null(width_in)) {
    width_in <- min(150, max(16, 0.06 * max(lengths(pages))))
  }

  starts  <- seq(1, n, by = page_size)
  ends    <- pmin(starts + page_size - 1, n)
  n_pages <- length(starts)

  ## --- determine device behaviour (current vs file) ---
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

  # flag to report overlap/clipping
  overlap_dropped <- FALSE

  ## --- inner page-drawing function ---
  draw_page <- function(x_from, x_to, page_index) {
    # For a single page, spread species across full width.
    # For multiple pages, use a fixed-width window per page so the last page is clumped.
    if (n_pages == 1L) {
      x_left  <- x_from
      x_right <- x_to
    } else {
      x_left  <- x_from
      x_right <- x_from + page_size - 1
    }

    ## add some horizontal padding so outer clusters are not on the border
    x_pad <- 0.5
    x_lim_left  <- x_left  - x_pad
    x_lim_right <- x_right + x_pad

    graphics::par(mar = c(8, 5, 1, 1), xaxs = "i", yaxs = "i")
    plot(c(x_lim_left, x_lim_right), c(0.1, -1),
         type = "n", xaxt = "n", yaxt = "n",
         xlab = "", ylab = expression(paste(phi, " coefficient")))
    graphics::axis(
      2, las = 2,
      at = seq(0, -1, by = -0.2),
      labels = seq(0, 1, by = 0.2)
    )

    # compute circle radius in data units once per page
    usr <- graphics::par("usr")  # c(x1, x2, y1, y2)
    pin <- graphics::par("pin")  # c(width_in_inch, height_in_inch)
    rx_d <- circle_inch * (usr[2] - usr[1]) / pin[1]
    ry_d <- circle_inch * (usr[4] - usr[3]) / pin[2]
    ymin <- usr[3]
    min_y <- ymin + ry_d  # centre so that bottom of circle touches axis

    # background bands
    if (!is.null(bands) && nrow(bands) > 0) {
      hit <- which(bands$x1 >= x_from & bands$x0 <= x_to)
      for (b in hit) {
        xl <- max(bands$x0[b], x_from) - 0.5
        xr <- min(bands$x1[b], x_to)   + 0.5
        graphics::rect(
          xl, bands$y0[b], xr, bands$y1[b],
          col = bands$col[b], border = bands$border[b]
        )
      }
    }

    # vertical legs to children
    for (r in 1:nrow(Pos)) {
      lx <- Pos[r, "lx"]; rx <- Pos[r, "rx"]
      if (lx >= x_from && lx <= x_to)
        graphics::segments(lx, Pos[r, "ly0"], lx, Pos[r, "ly1"])
      if (rx >= x_from && rx <= x_to)
        graphics::segments(rx, Pos[r, "ry0"], rx, Pos[r, "ry1"])
    }

    # horizontal connectors (merges)
    eps <- 1e-9
    hitH <- which(x1 >= x_from - eps & x0 <= x_to + eps)
    for (r in hitH) {
      hx0 <- max(x0[r], x_from)
      hx1 <- min(x1[r], x_to)
      if (hx1 >= hx0 - eps)
        graphics::segments(hx0, y[r], hx1, y[r])
    }

    # dashed cut line(s)
    if (!is.null(phi_cut)) {
      ycuts <- -phi_cut
      ycuts <- ycuts[is.finite(ycuts) & ycuts <= 0 & ycuts >= -1]
      if (length(ycuts)) {
        for (yy in ycuts) {
          graphics::segments(
            x_lim_left, yy, x_lim_right, yy,
            lty = 2, lwd = 1, col = "grey20"
          )
        }
      }
    }

    # tip labels (only for actual species on the page)
    idx <- x_from:x_to
    graphics::axis(
      1, at = idx,
      labels = species_names[species_order][idx],
      las = 2, cex.axis = cex_species
    )

    ## --- cluster / node labels ---
    if (identical(label_clusters, "all")) {

      ## 1) EXTRA: labels for phi = 0 using the same logic as "phi_cut" with phi_cut = 0
      ##    (we draw these FIRST and remember which cluster IDs we used)
      labels_phi0_drawn <- integer(0)

      if (!is.null(bands_phi0) && nrow(bands_phi0) > 0) {
        cross0 <- compute_phi_crossings(0, x_from, x_to, bands_phi0)
        if (!is.null(cross0) && nrow(cross0)) {
          # avoid clipping at bottom: lift centres to at least min_y
          low0 <- which(cross0$y < min_y)
          if (length(low0)) cross0$y[low0] <- min_y

          cross0_keep <- filter_non_overlapping(cross0, rx_d, ry_d)
          if (nrow(cross0_keep) < nrow(cross0)) {
            overlap_dropped <<- TRUE
          }

          if (nrow(cross0_keep) > 0) {
            # remember which cluster IDs we actually plotted at phi = 0
            labels_phi0_drawn <- unique(cross0_keep$label)

            graphics::symbols(
              cross0_keep$x, cross0_keep$y,
              circles = rep(1, nrow(cross0_keep)),
              inches  = circle_inch,
              add     = TRUE, bg = "white", fg = "black", lwd = 0.6
            )
            graphics::text(
              cross0_keep$x, cross0_keep$y,
              labels = cross0_keep$label,
              cex = cex_clusters
            )
          }
        }
      }

      ## 2) Standard "all nodes" labels:
      ##    where the parent branch meets the node (EXCLUDING clusters already
      ##    labelled at phi = 0 above)
      hit_nodes <- which(x1 >= x_from & x0 <= x_to)

      # do not re-label clusters that already got a phi = 0 label
      if (length(labels_phi0_drawn)) {
        hit_nodes <- setdiff(hit_nodes, labels_phi0_drawn)
      }

      if (length(hit_nodes)) {
        xm <- numeric(length(hit_nodes))
        ym <- numeric(length(hit_nodes))

        for (k in seq_along(hit_nodes)) {
          i_node <- hit_nodes[k]
          p_node <- parent_idx[i_node]

          if (!is.na(p_node)) {
            # node i_node is either left or right child of parent p_node
            if (CM[p_node, 1] == i_node) {
              x_lab <- Pos[p_node, "lx"]
            } else if (CM[p_node, 2] == i_node) {
              x_lab <- Pos[p_node, "rx"]
            } else {
              # safety fallback: use centre of this node
              x_lab <- center_x[i_node]
            }
            y_lab <- -H[p_node]  # parent height
          } else {
            # root: label at its own merge height, at centre of its span
            x_lab <- center_x[i_node]
            y_lab <- y[i_node]
          }

          xm[k] <- x_lab
          ym[k] <- y_lab
        }

        # keep only labels on current page horizontally
        keep <- which(xm >= x_from & xm <= x_to)
        if (length(keep)) {
          xm   <- xm[keep]
          ym   <- ym[keep]
          labs <- hit_nodes[keep]

          # avoid clipping at bottom: lift centres to at least min_y
          low <- which(ym < min_y)
          if (length(low)) ym[low] <- min_y

          df <- data.frame(x = xm, y = ym, label = labs, stringsAsFactors = FALSE)
          df_keep <- filter_non_overlapping(df, rx_d, ry_d)
          if (nrow(df_keep) < nrow(df)) {
            overlap_dropped <<- TRUE
          }

          if (nrow(df_keep) > 0) {
            graphics::symbols(
              df_keep$x, df_keep$y,
              circles = rep(1, nrow(df_keep)),
              inches  = circle_inch,
              add     = TRUE, bg = "white", fg = "black", lwd = 0.6
            )
            graphics::text(
              df_keep$x, df_keep$y,
              labels = df_keep$label,
              cex = cex_clusters
            )
          }
        }
      }

    } else if (identical(label_clusters, "phi_cut")) {
      if (!is.null(phi_cut) && !is.null(bands) && nrow(bands) > 0) {
        ycut <- -phi_cut
        for (yy in ycut) {
          cross <- compute_phi_crossings(yy, x_from, x_to, bands)
          if (!is.null(cross) && nrow(cross)) {
            # avoid clipping at bottom: lift centres to at least min_y
            low <- which(cross$y < min_y)
            if (length(low)) cross$y[low] <- min_y

            cross_keep <- filter_non_overlapping(cross, rx_d, ry_d)
            if (nrow(cross_keep) < nrow(cross)) {
              overlap_dropped <<- TRUE
            }

            if (nrow(cross_keep) > 0) {
              graphics::symbols(
                cross_keep$x, cross_keep$y,
                circles = rep(1, nrow(cross_keep)),
                inches  = circle_inch,
                add     = TRUE, bg = "white", fg = "black", lwd = 0.6
              )
              graphics::text(
                cross_keep$x, cross_keep$y,
                labels = cross_keep$label,
                cex = cex_clusters
              )
            }
          }
        }
      }
    }
  }

  ## --- current device: only first page ---
  if (use_current_dev) {
    p <- 1L
    draw_page(starts[p], ends[p], page_index = p)
    graphics::mtext(
      sprintf("Species %d-%d of %d", starts[p], ends[p], n),
      side = 3, line = 0.2, cex = 0.8
    )
    if (overlap_dropped) {
      warning(
        "Some cluster labels were omitted to avoid overlaps or clipping; ",
        "only non-overlapping cluster labels are plotted."
      )
    }
    return(invisible(NULL))
  }

  ## --- PDF output ---
  if (grepl("\\.pdf$", file, ignore.case = TRUE)) {
    grDevices::pdf(file, width = width_in, height = height_in, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)

    for (p in seq_along(starts)) {
      draw_page(starts[p], ends[p], page_index = p)
      graphics::mtext(
        sprintf("Species %d-%d of %d", starts[p], ends[p], n),
        side = 3, line = 0.2, cex = 0.8
      )
    }

    ## --- PNG output ---
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

      draw_page(starts[p], ends[p], page_index = p)
      graphics::mtext(
        sprintf("Species %d-%d of %d", starts[p], ends[p], n),
        side = 3, line = 0.2, cex = 0.8
      )

      grDevices::dev.off()
    }
  }

  if (overlap_dropped) {
    warning(
      "Some cluster labels were omitted to avoid overlaps or clipping; ",
      "only non-overlapping cluster labels are plotted."
    )
  }

  invisible(NULL)
}
