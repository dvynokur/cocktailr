#' Original Cocktail dendrogram (reference plot; minimal changes)
#'
#' Draw the Cocktail dendrogram using the plotting approach from Appendix S2
#' (Bruelheide 2016) with only minimal wrapper code. This mirrors the original:
#' species order via the "Species.sort" binary sequence and manual leg/crossbar segments.
#'
#' @param x A result list from [cocktail_cluster()] or [cocktail_cluster_orig()]
#'   containing `Cluster.species`, `Cluster.merged`, and `Cluster.height`.
#' @param highquality Integer: 1 uses a large PNG device (Cairo if available),
#'   0 draws to the current device. Default 1.
#' @param file Path to the output PNG when `highquality == 1`. Ignored otherwise.
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item `Species.sort` — character vector used to order species labels
#'   \item `Cluster.position` — matrix of segment coordinates
#' }
#'
#' @details
#' If **Cairo** is installed, the function uses the original `Cairo::Cairo()` call.
#' Otherwise it falls back to `grDevices::png()` with the same pixel dimensions.
#' The plotting logic itself is kept as in the original code block.
#'
#' @examples
#' \donttest{
#'   m <- as.matrix(vegan::dune); m[m > 0] <- 1L
#'   ref <- cocktailr::cocktail_cluster(m, progress = FALSE)
#'   plot_cocktail_orig_min(ref, highquality = 0)
#' }
#'
#' @seealso [cocktail_cluster()], [cocktail_cluster_orig()]
#' @export
#' @importFrom grDevices dev.off png
#' @importFrom graphics axis lines mtext par plot text

plot_cocktail_orig <- function(x, highquality = 1, file = "Cocktail_dendrogram.png") {
  stopifnot(is.list(x),
            !is.null(x$Cluster.species),
            !is.null(x$Cluster.merged),
            !is.null(x$Cluster.height))

  # ---- minimal glue: define the symbols the original code uses --------------
  Cluster.species <- x$Cluster.species
  Cluster.merged  <- x$Cluster.merged
  Cluster.height  <- x$Cluster.height

  # species names (the original code uses names(vegmatrix))
  sp_names <- colnames(Cluster.species)
  if (is.null(sp_names)) sp_names <- as.character(seq_len(ncol(Cluster.species)))

  # make a 1-row data.frame so that names(vegmatrix) works exactly like in original
  vegmatrix <- as.data.frame(matrix(0, nrow = 1, ncol = length(sp_names)))
  names(vegmatrix) <- sp_names

  n <- ncol(Cluster.species)  # number of species

  # ---- BEGIN: original plotting section (kept verbatim with only device guard) ----

  Species.sort <- array("", n)
  # holds the information in which sequence the species are arranged along the tips
  # of the Tree. Species.sort holds a string of binary numbers for all n-1 clusters
  # that show the assignment to a cluster in ascending order, i.e. from clusters formed last
  # to those formed first
  for (i in 1:n) {
    for (j in 1:n-1) {
      Species.sort[i] <- paste(Species.sort[i], Cluster.species[order(Cluster.height)[j], i], sep = "")
    }
  }

  highquality <- as.integer(highquality)
  # library(Cairo)  # original line; we avoid attaching, but use if available below

  if (highquality == 1) {
    if (requireNamespace("Cairo", quietly = TRUE)) {
      # exact original device if Cairo exists
      Cairo::Cairo(10976, 5664, file = file, type = "png", bg = "white")
    } else {
      # minimal fallback if Cairo is not installed (tiny deviation)
      grDevices::png(filename = file, width = 10976, height = 5664, bg = "white", res = NA)
    }
    # 10976, 5664
    #1372*8, 708*8, 8*8 standard size
    # 5488, 2832
    par(mar = c(50, 25, 1, 1)) # The default is c(5, 4, 4, 2) + 0.1.
    plot(c(1:(n)), c(seq(0.1, -1, (-0.1 - 1) / (n - 1))),
         type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.lab = 12)
    axis(1, las = 2, at = seq(1:n), labels = names(vegmatrix)[order(Species.sort)], cex.axis = 2)
    mtext(expression(paste(phi, " coefficient")), side = 2, line = 15, cex = 12)
    axis(2, las = 2, at = seq(0, -1, -0.2), labels = seq(0, 1, 0.2), cex.axis = 8)
  } else {
    plot(c(1:(n)), c(seq(0.1, -1, (-0.1 - 1) / (n - 1))),
         type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = expression(paste(phi, " coefficient")))
    axis(1, las = 2, at = seq(1:n), labels = names(vegmatrix)[order(Species.sort)], cex.axis = 1)
    axis(2, las = 2, at = seq(0, -1, -0.2), labels = seq(0, 1, 0.2))
  }

  Cluster.position <- array(
    NA, c(n - 1, 6),
    dimnames = list(c(1:(n - 1)),
                    c("left.leg.x", "left.leg.y0", "left.leg.y1",
                      "right.leg.x", "right.leg.y0", "right.leg.y1"))
  )
  # Array with x and y coordinates for all clusters. Coordinates have a negative
  # sign to arrange the tree in descending values of the phi.index
  # y0: lower part of the leg
  # y1: upper part of the leg
  for (i in 1:(n - 1)) {
    # loop for all clusters
    if (Cluster.merged[i, 1] < 0) {
      # left leg
      Cluster.position[i, 1] <- which(order(Species.sort) == -Cluster.merged[i, 1])
      Cluster.position[i, 2] <- -1
      Cluster.position[i, 3] <- -Cluster.height[i]
    } else {
      Cluster.position[i, 1] <-
        mean(match(which(Cluster.species[Cluster.merged[i, 1], ] == 1), order(Species.sort)))
      Cluster.position[i, 2] <- -Cluster.height[Cluster.merged[i, 1]]
      Cluster.position[i, 3] <- -Cluster.height[i]
    }
    if (Cluster.merged[i, 2] < 0) {
      # right leg
      Cluster.position[i, 4] <- which(order(Species.sort) == -Cluster.merged[i, 2])
      Cluster.position[i, 5] <- -1
      Cluster.position[i, 6] <- -Cluster.height[i]
    } else {
      Cluster.position[i, 4] <-
        mean(match(which(Cluster.species[Cluster.merged[i, 2], ] == 1), order(Species.sort)))
      Cluster.position[i, 5] <- -Cluster.height[Cluster.merged[i, 2]]
      Cluster.position[i, 6] <- -Cluster.height[i]
    }
    # labels
    text(mean(c(Cluster.position[i, 1], Cluster.position[i, 4])),
         Cluster.position[i, 6] - 0.01, i, cex = 1)
    # left leg
    lines(c(Cluster.position[i, 1], Cluster.position[i, 1]),
          c(Cluster.position[i, 2], Cluster.position[i, 3]))
    # right leg
    lines(c(Cluster.position[i, 4], Cluster.position[i, 4]),
          c(Cluster.position[i, 5], Cluster.position[i, 6]))
    # horizontal line
    lines(c(Cluster.position[i, 1], Cluster.position[i, 4]),
          c(Cluster.position[i, 6], Cluster.position[i, 6]))
  }
  if (highquality == 1) {
    grDevices::dev.off() # creates a file with the above plot
    par(mar = c(5, 4, 4, 2)) # The default is c(5, 4, 4, 2) + 0.1.
  }

  # ---- END original plotting section ---------------------------------------

  invisible(list(Species.sort = Species.sort, Cluster.position = Cluster.position))
}
