#' Plot a Cocktail dendrogram (PDF/PNG or current device)
#'
#' @description
#' Draw a dendrogram-style plot of a fitted Cocktail clustering.
#' Species (tips) are arranged left-to-right in the same order used
#' by the Cocktail algorithm, and the y-axis shows the \eqn{\phi}
#' height of each merge.
#'
#' Optionally:
#' * highlight “parent” clusters at or above a chosen \eqn{\phi} cut
#'   with coloured bands, and
#' * label internal nodes or only the parent clusters at the cut.
#'
#' The plot can be:
#' * written to a **PDF** (multi-page) when `file` ends with `.pdf`,
#' * written to one or several **PNG** files when `file` ends with `.png`, or
#' * drawn on the current graphics device (e.g. RStudio Plots pane) when
#'   `file` is `NULL` or has no recognized extension. In this “interactive”
#'   mode, only the **first page** is drawn.
#'
#' @param x A Cocktail object (e.g. from [cocktail_cluster()]) with
#'   components `Cluster.species`, `Cluster.merged`, and `Cluster.height`.
#'   If `x$species` is present it is used for tip labels; otherwise
#'   `colnames(x$Cluster.species)` are used.
#' @param file Optional path to output file.
#'   * If it ends with `.pdf`, a multi-page PDF is created.
#'   * If it ends with `.png`, one or several PNG files are created.
#'   * If `NULL` or the extension is not `.pdf`/`.png`, no file is written
#'     and the plot is drawn on the current device (first page only).
#' @param phi_cut Optional numeric in (0,1). If given, parent clusters
#'   at or above this \eqn{\phi} are highlighted with coloured bands and
#'   a dashed horizontal cut line at `phi_cut` is drawn.
#' @param label_clusters Logical/character. One of:
#'   * `FALSE` (default): no node labels;
#'   * `"all"`: label **every** internal node (merge);
#'   * `"phi_cut"`: label clusters that intersect the `phi_cut` line
#'     (requires `phi_cut` to be non-`NULL`).
#'   Labels are node indices (row numbers of `Cluster.species`), shown
#'   without any `"c_"` prefix.
#' @param cex_species Numeric expansion factor for **species (tip) labels**
#'   on the x-axis.
#' @param cex_clusters Numeric expansion factor for **cluster-number labels**
#'   (node bubbles).
#' @param circle_inch Radius (inches) of the white/black bubble drawn behind
#'   cluster labels.
#' @param page_size Integer, number of tips per page. Use
#'   `ncol(x$Cluster.species)` to place all tips on a single page.
#' @param width_in Plot width in inches. If `NULL`, a width is chosen
#'   automatically (capped at 150 inches) based on `page_size`.
#' @param height_in Plot height in inches (default `10`).
#' @param alpha_fill,alpha_border Numeric alpha levels (0–1) for band fill
#'   and border, respectively.
#' @param palette Color palette name. If it matches one of
#'   `grDevices::hcl.pals()`, colours are taken from
#'   `grDevices::hcl.colors()`. If `"rainbow"` (default, case-insensitive),
#'   colours are taken from `grDevices::rainbow()`. Otherwise the
#'   fallback palette `"Set3"` is used.
#' @param png_res Resolution in dpi for PNG output (ignored for PDF).
#'
#' @details
#' Species order is reproduced from the Cocktail implementation by
#' sorting the per-tip membership strings derived from `x$Cluster.species`.
#'
#' Large trees are automatically paginated according to `page_size`. When
#' writing to PDF, all pages go into a single file. When writing to PNG,
#' multiple pages are written as `base_p01.png`, `base_p02.png`, etc.,
#' where `base` is the filename without the `.png` extension.
#'
#' When `file` is `NULL` (or has an unrecognized extension), only the
#' first page is drawn on the **current** graphics device, which is
#' typically the Plots pane in RStudio.
#'
#' @return
#' Invisibly returns `NULL`. The side effect is the creation of a PDF/PNG
#' on disk (if `file` is a supported filename) and/or drawing the plot
#' on the active device.
#'
#' @importFrom graphics axis par segments rect mtext plot symbols text
#' @importFrom grDevices pdf png rainbow hcl.colors hcl.pals adjustcolor
#' @importFrom tools file_ext
#' @export
cocktail_plot <- function(
    x,
    file           = NULL,
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
  ## --- basic checks --------------------------------------------------------
  if (!is.list(x) ||
      !all(c("Cluster.species", "Cluster.merged", "Cluster.height") %in% names(x))) {
    stop("`x` must be a Cocktail object with components: ",
         "Cluster.species, Cluster.merged, Cluster.height.")
  }

  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height

  if (!is.matrix(CS) || !is.matrix(CM) || !is.numeric(H)) {
    stop("`Cluster.species` and `Cluster.merged` must be matrices, ",
         "`Cluster.height` must be numeric.")
  }

  n <- ncol(CS)
  if (n < 2L) stop("Need at least 2 species (columns in Cluster.species).")

  # Species names: prefer x$species, else colnames(CS)
  species_names <- if (!is.null(x$species)) {
    x$species
  } else {
    colnames(CS)
  }
  if (is.null(species_names) || length(species_names) != n) {
    species_names <- paste0("sp_", seq_len(n))
  }

  ## --- species order (Cocktail-style) --------------------------------------
  ord <- order(H)
  Species.sort <- vapply(
    X   = seq_len(n),
    FUN = function(i) paste(CS[ord, i], collapse = ""),
    FUN.VALUE = character(1)
  )
  species_order <- order(Species.sort)

  ## --- precompute coordinates for every merge ------------------------------
  Pos <- matrix(
    NA_real_,
    nrow = n - 1L, ncol = 6L,
    dimnames = list(
      as.character(seq_len(n - 1L)),
      c("lx", "ly0", "ly1", "rx", "ry0", "ry1")
    )
  )

  for (i in seq_len(n - 1L)) {
    ## left child
    if (CM[i, 1L] < 0L) {
      tip <- -CM[i, 1L]
      Pos[i, "lx"]  <- which(species_order == tip)
      Pos[i, "ly0"] <- -1
      Pos[i, "ly1"] <- -H[i]
    } else {
      spp <- which(CS[CM[i, 1L], ] == 1L)
      Pos[i, "lx"]  <- mean(match(spp, species_order))
      Pos[i, "ly0"] <- -H[CM[i, 1L]]
      Pos[i, "ly1"] <- -H[i]
    }

    ## right child
    if (CM[i, 2L] < 0L) {
      tip <- -CM[i, 2L]
      Pos[i, "rx"]  <- which(species_order == tip)
      Pos[i, "ry0"] <- -1
      Pos[i, "ry1"] <- -H[i]
    } else {
      spp <- which(CS[CM[i, 2L], ] == 1L)
      Pos[i, "rx"]  <- mean(match(spp, species_order))
      Pos[i, "ry0"] <- -H[CM[i, 2L]]
      Pos[i, "ry1"] <- -H[i]
    }
  }

  x0 <- pmin(Pos[, "lx"], Pos[, "rx"])
  x1 <- pmax(Pos[, "lx"], Pos[, "rx"])
  y  <- -H[seq_len(n - 1L)]

  ## --- bands for parent clusters at phi_cut --------------------------------
  bands   <- NULL
  top_idx <- integer(0L)

  if (!is.null(phi_cut)) {
    idx <- which(H >= phi_cut)
    if (length(idx)) {
      children <- unique(as.integer(CM[idx, , drop = FALSE][CM[idx, , drop = FALSE] > 0]))
      top_idx  <- sort(setdiff(idx, intersect(idx, children)))

      if (length(top_idx)) {
        k <- length(top_idx)

        cols_base <- if (palette %in% rownames(grDevices::hcl.pals())) {
          grDevices::hcl.colors(max(k, 3L), palette = palette)
        } else if (tolower(palette) == "rainbow") {
          grDevices::rainbow(k)
        } else {
          grDevices::hcl.colors(max(k, 3L), palette = "Set3")
        }
        cols_base <- cols_base[seq_len(k)]

        blist <- vector("list", k)
        for (ii in seq_along(top_idx)) {
          i <- top_idx[ii]
          spp  <- which(CS[i, ] == 1L)
          tips <- sort(match(spp, species_order))
          if (length(tips) < 2L) next
          blist[[ii]] <- data.frame(
            x0 = min(tips),
            x1 = max(tips),
            y0 = -1,
            y1 = 0.1,
            col    = grDevices::adjustcolor(cols_base[ii], alpha.f = alpha_fill),
            border = grDevices::adjustcolor(cols_base[ii], alpha.f = alpha_border),
            node_id = i,
            stringsAsFactors = FALSE
          )
        }
        bands <- do.call(rbind, blist)
      }
    }
  }

  ## --- helper: compute cluster labels at phi_cut ---------------------------
  compute_phi_crossings <- function(ycut, x_from, x_to) {
    if (is.null(bands) || nrow(bands) == 0L) return(NULL)
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
        bcand <- which(bands$x0 <= x_leg & bands$x1 >= x_leg)
        if (!length(bcand)) next
        nid <- bands$node_id[bcand[1L]]
        crosses[[length(crosses) + 1L]] <- data.frame(
          x = x_leg,
          y = ycut,
          label = nid,
          stringsAsFactors = FALSE
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
        bcand <- which(bands$x0 <= x_leg & bands$x1 >= x_leg)
        if (!length(bcand)) next
        nid <- bands$node_id[bcand[1L]]
        crosses[[length(crosses) + 1L]] <- data.frame(
          x = x_leg,
          y = ycut,
          label = nid,
          stringsAsFactors = FALSE
        )
      }
    }

    if (!length(crosses)) return(NULL)
    do.call(rbind, crosses)
  }

  ## --- pagination and width ------------------------------------------------
  pages <- split(seq_len(n), ceiling(seq_len(n) / page_size))

  if (is.null(width_in)) {
    width_in <- min(150, max(16, 0.06 * max(lengths(pages))))
  }

  ## --- inner page drawing function ----------------------------------------
  draw_page <- function(x_from, x_to) {
    graphics::par(mar = c(8, 5, 1, 1), xaxs = "i", yaxs = "i")

    graphics::plot(
      c(x_from, x_to), c(0.1, -1),
      type = "n", xaxt = "n", yaxt = "n",
      xlab = "", ylab = expression(paste(phi, " coefficient"))
    )
    graphics::axis(
      side = 2, las = 2,
      at = seq(0, -1, by = -0.2),
      labels = seq(0, 1, by = 0.2)
    )

    ## background bands (clip to page; +/-0.5 to cover tick centres)
    if (!is.null(bands) && nrow(bands) > 0L) {
      hit <- which(bands$x1 >= x_from & bands$x0 <= x_to)
      for (b in hit) {
        xl <- max(bands$x0[b], x_from) - 0.5
        xr <- min(bands$x1[b], x_to)   + 0.5
        graphics::rect(
          xl, bands$y0[b], xr, bands$y1[b],
          col = bands$col[b],
          border = bands$border[b]
        )
      }
    }

    ## vertical legs
    for (r in seq_len(nrow(Pos))) {
      lx <- Pos[r, "lx"]; rx <- Pos[r, "rx"]
      if (lx >= x_from && lx <= x_to) {
        graphics::segments(lx, Pos[r, "ly0"], lx, Pos[r, "ly1"])
      }
      if (rx >= x_from && rx <= x_to) {
        graphics::segments(rx, Pos[r, "ry0"], rx, Pos[r, "ry1"])
      }
    }

    ## horizontal connectors (robust clipping)
    eps  <- 1e-9
    hitH <- which(x1 >= x_from - eps & x0 <= x_to + eps)
    for (r in hitH) {
      hx0 <- max(x0[r], x_from)
      hx1 <- min(x1[r], x_to)
      if (hx1 >= hx0 - eps) {
        graphics::segments(hx0, y[r], hx1, y[r])
      }
    }

    ## dashed cut line at phi_cut
    if (!is.null(phi_cut)) {
      yy <- -phi_cut
      if (is.finite(yy) && yy <= 0 && yy >= -1) {
        graphics::segments(
          x_from - 0.5, yy,
          x_to   + 0.5, yy,
          lty = 2, lwd = 1, col = "grey20"
        )
      }
    }

    ## tip labels
    idx <- x_from:x_to
    graphics::axis(
      side = 1,
      at   = idx,
      labels = species_names[species_order][idx],
      las = 2,
      cex.axis = cex_species
    )

    ## cluster labels --------------------------------------------------------
    if (identical(label_clusters, "all")) {
      hit_nodes <- which(x1 >= x_from & x0 <= x_to)
      if (length(hit_nodes)) {
        xm <- (pmax(x0[hit_nodes], x_from) + pmin(x1[hit_nodes], x_to)) / 2
        ym <- y[hit_nodes]

        graphics::symbols(
          xm, ym,
          circles = rep(1, length(xm)),
          inches  = circle_inch,
          add     = TRUE,
          bg = "white", fg = "black", lwd = 0.6
        )
        graphics::text(
          xm, ym,
          labels = hit_nodes,
          cex    = cex_clusters
        )
      }

    } else if (identical(label_clusters, "phi_cut") && !is.null(phi_cut) && length(top_idx)) {
      yy <- -phi_cut
      cross <- compute_phi_crossings(yy, x_from, x_to)
      if (!is.null(cross) && nrow(cross)) {
        graphics::symbols(
          cross$x, cross$y,
          circles = rep(1, nrow(cross)),
          inches  = circle_inch,
          add     = TRUE,
          bg = "white", fg = "black", lwd = 0.6
        )
        graphics::text(
          cross$x, cross$y,
          labels = cross$label,
          cex    = cex_clusters
        )
      }
    }
  }

  starts <- seq(1L, n, by = page_size)
  ends   <- pmin(starts + page_size - 1L, n)

  ## --- device handling -----------------------------------------------------
  if (is.null(file)) {
    # interactive mode: draw only the first page on current device
    draw_page(starts[1L], ends[1L])
    graphics::mtext(
      sprintf("Species %d-%d of %d", starts[1L], ends[1L], n),
      side = 3, line = 0.2, cex = 0.8
    )
    return(invisible(NULL))
  }

  ext <- tolower(tools::file_ext(file))

  if (!ext %in% c("pdf", "png")) {
    warning(
      "Unrecognized file extension for `file` (expected .pdf or .png). ",
      "No file will be written; drawing first page on the current device.\n",
      "To save the plot, please use a filename ending with '.pdf' or '.png'."
    )
    draw_page(starts[1L], ends[1L])
    graphics::mtext(
      sprintf("Species %d-%d of %d", starts[1L], ends[1L], n),
      side = 3, line = 0.2, cex = 0.8
    )
    return(invisible(NULL))
  }

  if (ext == "pdf") {
    grDevices::pdf(file, width = width_in, height = height_in, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)

    for (p in seq_along(starts)) {
      draw_page(starts[p], ends[p])
      graphics::mtext(
        sprintf("Species %d-%d of %d", starts[p], ends[p], n),
        side = 3, line = 0.2, cex = 0.8
      )
    }

  } else if (ext == "png") {
    base <- sub("\\.png$", "", file, ignore.case = TRUE)
    n_pages <- length(starts)

    for (p in seq_along(starts)) {
      fname <- if (n_pages > 1L) {
        sprintf("%s_p%02d.png", base, p)
      } else {
        file
      }
      grDevices::png(
        filename = fname,
        width  = width_in,
        height = height_in,
        units  = "in",
        res    = png_res
      )
      draw_page(starts[p], ends[p])
      graphics::mtext(
        sprintf("Species %d-%d of %d", starts[p], ends[p], n),
        side = 3, line = 0.2, cex = 0.8
      )
      grDevices::dev.off()
    }
  }

  invisible(NULL)
}
