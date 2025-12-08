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
#' @param phi_cut Optional numeric. If given (and \code{clusters} is
#'   \code{NULL}), draw background bands for parent clusters with
#'   \eqn{\phi \ge \mathrm{phi\_cut}}, and a dashed horizontal line
#'   at this \eqn{\phi} level across each page.
#' @param label_clusters Logical. If \code{FALSE}, no internal node labels
#'   are drawn. If \code{TRUE}:
#'   \itemize{
#'     \item if \code{clusters} is supplied, label those cluster nodes
#'           (topmost per group) at their elbows, except that any supplied
#'           clusters whose vertical branch crosses \code{y = 0} are labelled
#'           \emph{only once} at \code{y = 0};
#'     \item else if \code{phi_cut} is supplied, label the clusters at
#'           the \eqn{\phi\_cut} level (where vertical branches cross the
#'           dashed line);
#'     \item else, label all internal nodes.
#'   }
#'   Labels are numeric node IDs.
#' @param clusters Optional selection of clusters (nodes) to be shown as
#'   coloured bands and used for labelling when \code{label_clusters = TRUE}.
#'   Can be:
#'   \itemize{
#'     \item a vector like \code{c("c_12","c_27")} or \code{c(12,27)},
#'           where each element defines a separate group; or
#'     \item a \strong{list} of such vectors, where each element defines
#'           a \strong{union group} of nodes (e.g.
#'           \code{list(c(1,4,6,14), c(2,12,17))} or your
#'           \code{node_groups <- split(names(grp_nodes), grp_nodes)}).
#'   }
#'   Within each group, if both an ancestor and descendant are present,
#'   only the topmost (ancestor) is used for elbow labelling and band
#'   placement. All valid supplied IDs are considered for labels at
#'   \eqn{\phi = 0} when their branch crosses that level.
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
    clusters       = NULL,
    cex_species    = 0.3,
    cex_clusters   = 0.6,
    circle_inch    = 0.1,
    page_size      = 300,
    width_in       = NULL,
    height_in      = 10,
    alpha_fill     = 0.18,
    alpha_border   = 0.65,
    palette        = "rainbow",
    png_res        = 300,
    main           = NULL,
    ...
) {
  ## ---- helpers ------------------------------------------------------------
  .parse_clusters_arg <- function(v) {
    if (is.null(v)) return(NULL)

    parse_one <- function(x) {
      if (is.character(x)) {
        as.integer(sub("^c_", "", x))
      } else {
        as.integer(x)
      }
    }

    if (is.list(v)) {
      lapply(v, parse_one)
    } else {
      ids <- parse_one(v)
      lapply(as.list(ids), identity)
    }
  }

  .keep_topmost_within <- function(ids, CM) {
    ids <- sort(unique(ids))
    if (!length(ids)) return(ids)
    cm_sub <- CM[ids, , drop = FALSE]
    kids <- unique(as.integer(cm_sub[cm_sub > 0]))
    sort(setdiff(ids, intersect(ids, kids)))
  }

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

  compute_phi_crossings <- function(ycut, x_from, x_to, bands_df, Pos) {
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
        nid <- bands_df$node_id[bcand[1L]]
        crosses[[length(crosses) + 1L]] <- data.frame(
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
        nid <- bands_df$node_id[bcand[1L]]
        crosses[[length(crosses) + 1L]] <- data.frame(
          x = x_leg, y = ycut, label = nid
        )
      }
    }

    if (!length(crosses)) return(NULL)
    do.call(rbind, crosses)
  }

  ## small helper to choose page title
  .page_title <- function(p, starts, ends, n, main) {
    if (is.null(main)) {
      sprintf("Species %d-%d of %d", starts[p], ends[p], n)
    } else {
      if (length(main) >= p) {
        main[p]
      } else {
        main[1L]
      }
    }
  }

  ## ---- basic setup --------------------------------------------------------
  if (!is.null(clusters) && !is.null(phi_cut)) {
    warning("Both `clusters` and `phi_cut` supplied; using `clusters` for bands and labels.")
  }

  CS <- x$Cluster.species
  CM <- x$Cluster.merged
  H  <- x$Cluster.height
  n  <- ncol(CS)

  species_names <- if (!is.null(x$species)) x$species else colnames(CS)

  ## ---- species order (Cocktail-like) --------------------------------------
  ord <- order(H)
  Species.sort <- vapply(
    seq_len(n),
    function(i) paste(CS[ord, i], collapse = ""),
    FUN.VALUE = character(1)
  )
  species_order <- order(Species.sort)

  ## ---- precompute coordinates for every merge (legs + y) ------------------
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

  ## ---- cluster "centre" (horizontal midpoint of tips) ---------------------
  center_x <- (x0 + x1) / 2

  ## ---- parent index for each node (1..n-1), NA for roots ------------------
  parent_idx <- rep(NA_integer_, nrow(CM))
  for (p in seq_len(nrow(CM))) {
    for (j in 1:2) {
      ch <- CM[p, j]
      if (ch > 0L && ch <= nrow(CM)) {
        parent_idx[ch] <- p
      }
    }
  }

  ## For each node: which parent row/leg (to know its vertical branch)
  child_parent_row <- parent_idx
  child_leg_side   <- rep(NA_character_, length(parent_idx))
  for (p in seq_len(nrow(CM))) {
    if (CM[p, 1] > 0L) child_leg_side[CM[p, 1]] <- "L"
    if (CM[p, 2] > 0L) child_leg_side[CM[p, 2]] <- "R"
  }

  ## ---- bands: phi_cut-based (if no clusters) ------------------------------
  bands_phi  <- NULL

  if (is.null(clusters) && !is.null(phi_cut)) {
    idx <- which(H >= phi_cut)
    if (length(idx)) {
      children <- unique(as.integer(CM[idx, , drop = FALSE][CM[idx, , drop = FALSE] > 0]))
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
        bands_phi <- do.call(rbind, blist)
      }
    }
  }

  ## ---- bands: custom cluster groups (if clusters supplied) ----------------
  clusters_top_nodes   <- NULL
  bands_custom         <- NULL
  clusters_all_ids_val <- NULL  # all valid cluster IDs supplied (for y=0 labels)

  if (!is.null(clusters)) {
    n_nodes <- nrow(CS)
    grp_ids <- .parse_clusters_arg(clusters)

    if (length(grp_ids)) {
      # validate and drop invalid IDs per group
      valid <- logical(length(grp_ids))
      for (g in seq_along(grp_ids)) {
        ids <- grp_ids[[g]]
        bad <- is.na(ids) | ids < 1L | ids > n_nodes
        if (any(bad)) {
          warning(
            "Dropping invalid node IDs in custom cluster group ", g, ": ",
            paste(ids[bad], collapse = ", ")
          )
          ids <- ids[!bad]
        }
        grp_ids[[g]] <- ids
        valid[g] <- length(ids) > 0L
      }
      grp_ids <- grp_ids[valid]

      if (length(grp_ids)) {
        # all valid IDs (for y=0 labels)
        clusters_all_ids_val <- sort(unique(unlist(grp_ids)))

        # keep only topmost nodes per group (for bands + elbow labels)
        top_nodes <- lapply(grp_ids, .keep_topmost_within, CM = CM)

        # build bands over union of species of topmost nodes
        bx0 <- numeric(0)
        bx1 <- numeric(0)

        for (g in seq_along(top_nodes)) {
          ids_top <- top_nodes[[g]]
          if (!length(ids_top)) next
          sp_logical <- colSums(CS[ids_top, , drop = FALSE] > 0) > 0L
          spp <- which(sp_logical)
          if (length(spp) < 2) next
          tips <- sort(match(spp, species_order))
          bx0 <- c(bx0, min(tips))
          bx1 <- c(bx1, max(tips))
        }

        if (length(bx0)) {
          bands_custom <- data.frame(
            x0 = bx0,
            x1 = bx1,
            y0 = -1,
            y1 = 0.1,
            stringsAsFactors = FALSE
          )

          k <- nrow(bands_custom)
          cols_base <-
            if (palette %in% rownames(grDevices::hcl.pals())) {
              grDevices::hcl.colors(max(k, 3), palette = palette)
            } else if (tolower(palette) == "rainbow") {
              grDevices::rainbow(k)
            } else {
              grDevices::hcl.colors(max(k, 3), palette = "Set3")
            }
          cols_base <- cols_base[seq_len(k)]

          bands_custom$col    <- grDevices::adjustcolor(cols_base, alpha.f = alpha_fill)
          bands_custom$border <- grDevices::adjustcolor(cols_base, alpha.f = alpha_border)

          clusters_top_nodes <- top_nodes
        }
      }
    }
  }

  ## ---- pagination (+ auto width if needed) --------------------------------
  pages <- split(seq_len(n), ceiling(seq_len(n) / page_size))
  if (is.null(width_in)) {
    width_in <- min(150, max(16, 0.06 * max(lengths(pages))))
  }

  starts  <- seq(1, n, by = page_size)
  ends    <- pmin(starts + page_size - 1, n)
  n_pages <- length(starts)

  ## ---- determine device behaviour (current vs file) -----------------------
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

  ## ---- inner page-drawing function ----------------------------------------
  draw_page <- function(x_from, x_to, page_index, ...) {
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
         xlab = "", ylab = expression(paste(phi, " coefficient")),
         ...)

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

    # choose which bands to draw (custom takes precedence)
    bands_to_draw <- NULL
    if (!is.null(bands_custom) && nrow(bands_custom) > 0) {
      bands_to_draw <- bands_custom
    } else if (!is.null(bands_phi) && nrow(bands_phi) > 0) {
      bands_to_draw <- bands_phi
    }

    # background bands
    if (!is.null(bands_to_draw) && nrow(bands_to_draw) > 0) {
      hit <- which(bands_to_draw$x1 >= x_from & bands_to_draw$x0 <= x_to)
      for (b in hit) {
        xl <- max(bands_to_draw$x0[b], x_from) - 0.5
        xr <- min(bands_to_draw$x1[b], x_to)   + 0.5
        graphics::rect(
          xl, bands_to_draw$y0[b], xr, bands_to_draw$y1[b],
          col    = bands_to_draw$col[b],
          border = bands_to_draw$border[b]
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

    ## ---- cluster / node labels -------------------------------------------

    if (!isTRUE(label_clusters)) return(invisible(NULL))

    ## Case A: clusters supplied -> labels for those clusters
    if (!is.null(clusters_top_nodes) && length(clusters_top_nodes)) {

      ## A1. First: labels at y = 0 for supplied clusters whose branch crosses y = 0
      ids_crossing0_page <- integer(0)
      df0_keep <- NULL

      if (!is.null(clusters_all_ids_val) && length(clusters_all_ids_val)) {
        x_phi0 <- numeric(0)
        y_phi0 <- numeric(0)
        lab0   <- integer(0)

        for (i_node in clusters_all_ids_val) {
          if (is.na(i_node) || i_node < 1L || i_node > length(child_parent_row))
            next
          p <- child_parent_row[i_node]
          if (is.na(p)) next  # no parent -> no vertical branch

          side <- child_leg_side[i_node]
          if (side == "L") {
            x_leg  <- Pos[p, "lx"]
            y0_leg <- Pos[p, "ly0"]
            y1_leg <- Pos[p, "ly1"]
          } else if (side == "R") {
            x_leg  <- Pos[p, "rx"]
            y0_leg <- Pos[p, "ry0"]
            y1_leg <- Pos[p, "ry1"]
          } else {
            next
          }

          # does this vertical branch cross y=0?
          if (min(y0_leg, y1_leg) <= 0 && max(y0_leg, y1_leg) >= 0) {
            if (x_leg < x_from || x_leg > x_to) next
            x_phi0 <- c(x_phi0, x_leg)
            y_phi0 <- c(y_phi0, 0)
            lab0   <- c(lab0, i_node)
          }
        }

        if (length(x_phi0)) {
          df0 <- data.frame(
            x = x_phi0,
            y = y_phi0,
            label = lab0,
            stringsAsFactors = FALSE
          )

          df0_keep <- filter_non_overlapping(df0, rx_d, ry_d)
          if (nrow(df0_keep) < nrow(df0)) {
            overlap_dropped <<- TRUE
          }

          if (nrow(df0_keep) > 0) {
            ids_crossing0_page <- unique(df0_keep$label)

            graphics::symbols(
              df0_keep$x, df0_keep$y,
              circles = rep(1, nrow(df0_keep)),
              inches  = circle_inch,
              add     = TRUE, bg = "white", fg = "black", lwd = 0.6
            )
            graphics::text(
              df0_keep$x, df0_keep$y,
              labels = df0_keep$label,
              cex = cex_clusters
            )
          }
        }
      }

      ## A2. Elbow labels for topmost nodes, but NOT for those already labelled at y=0
      node_ids <- sort(unique(unlist(clusters_top_nodes)))
      node_ids <- node_ids[node_ids >= 1L & node_ids <= length(parent_idx)]
      if (length(ids_crossing0_page)) {
        node_ids <- setdiff(node_ids, ids_crossing0_page)
      }

      if (length(node_ids)) {
        xm <- numeric(length(node_ids))
        ym <- numeric(length(node_ids))

        for (k in seq_along(node_ids)) {
          i_node <- node_ids[k]
          p_node <- parent_idx[i_node]

          if (!is.na(p_node)) {
            # label where this cluster begins: parent horizontal -> vertical
            if (CM[p_node, 1] == i_node) {
              x_lab <- Pos[p_node, "lx"]
            } else if (CM[p_node, 2] == i_node) {
              x_lab <- Pos[p_node, "rx"]
            } else {
              x_lab <- center_x[i_node]
            }
            y_lab <- -H[p_node]
          } else {
            # root
            x_lab <- center_x[i_node]
            y_lab <- y[i_node]
          }

          xm[k] <- x_lab
          ym[k] <- y_lab
        }

        keep <- which(xm >= x_from & xm <= x_to & !is.na(xm) & !is.na(ym))
        if (length(keep)) {
          df <- data.frame(
            x = xm[keep],
            y = ym[keep],
            label = node_ids[keep],
            stringsAsFactors = FALSE
          )

          low <- which(df$y < min_y)
          if (length(low)) df$y[low] <- min_y

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

      return(invisible(NULL))
    }

    ## Case B: no clusters, but phi_cut present -> label φ-cut clusters
    if (!is.null(phi_cut) && !is.null(bands_phi) && nrow(bands_phi) > 0) {
      ycut <- -phi_cut
      cross <- compute_phi_crossings(ycut, x_from, x_to, bands_phi, Pos)
      if (!is.null(cross) && nrow(cross)) {
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

      return(invisible(NULL))
    }

    ## Case C: no clusters, no phi_cut -> label all internal nodes
    hit_nodes <- which(x1 >= x_from & x0 <= x_to)
    if (!length(hit_nodes)) return(invisible(NULL))

    xm <- numeric(length(hit_nodes))
    ym <- numeric(length(hit_nodes))

    for (k in seq_along(hit_nodes)) {
      i_node <- hit_nodes[k]
      p_node <- parent_idx[i_node]

      if (!is.na(p_node)) {
        if (CM[p_node, 1] == i_node) {
          x_lab <- Pos[p_node, "lx"]
        } else if (CM[p_node, 2] == i_node) {
          x_lab <- Pos[p_node, "rx"]
        } else {
          x_lab <- center_x[i_node]
        }
        y_lab <- -H[p_node]
      } else {
        x_lab <- center_x[i_node]
        y_lab <- y[i_node]
      }

      xm[k] <- x_lab
      ym[k] <- y_lab
    }

    keep <- which(xm >= x_from & xm <= x_to)
    if (!length(keep)) return(invisible(NULL))

    df <- data.frame(
      x = xm[keep],
      y = ym[keep],
      label = hit_nodes[keep],
      stringsAsFactors = FALSE
    )

    low <- which(df$y < min_y)
    if (length(low)) df$y[low] <- min_y

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

    invisible(NULL)
  }

  ## ---- current device: only first page ------------------------------------
  if (use_current_dev) {
    p <- 1L
    draw_page(starts[p], ends[p], page_index = p, ...)
    graphics::mtext(
      .page_title(p, starts, ends, n, main),
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

  ## ---- PDF output ---------------------------------------------------------
  if (grepl("\\.pdf$", file, ignore.case = TRUE)) {
    grDevices::pdf(file, width = width_in, height = height_in, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)

    for (p in seq_along(starts)) {
      draw_page(starts[p], ends[p], page_index = p, ...)
      graphics::mtext(
        .page_title(p, starts, ends, n, main),
        side = 3, line = 0.2, cex = 0.8
      )
    }

    ## ---- PNG output ---------------------------------------------------------
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

      draw_page(starts[p], ends[p], page_index = p, ...)
      graphics::mtext(
        .page_title(p, starts, ends, n, main),
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
