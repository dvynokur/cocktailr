#' Dendrogram of fuzzy Cocktail groups (groups as tips)
#'
#' @param fuzzy  data.frame from cocktail_fuzzy(); columns must include:
#'               group, species, phi, cut
#' @param cut    numeric phi cut to use (e.g. 0.30)
#' @param min_phi keep only memberships with phi >= min_phi (default 0)
#' @param dist   distance between groups: "bray" (default) or "cosine"
#' @param hclust_method linkage for hclust: "average" (UPGMA), "complete", etc.
#' @param group_names optional named character vector mapping group codes to nice names
#'                    e.g. c(G107="Stipa lessingiana grp", G114="Filipendula grp")
#' @param pdf_file optional path; if given, plot is written to this PDF
#' @param width_in,height_in PDF size if pdf_file is used (width auto-scales if NULL)
#' @return list(hclust=, dist=, X=matrix, labels=character)
#' @export
fuzzy_groups_dendrogram <- function(
    fuzzy,
    cut,
    min_phi = 0,
    dist = c("bray","cosine"),
    hclust_method = "average",
    group_names = NULL,
    pdf_file = NULL,
    width_in = NULL,
    height_in = 7
){
  dist <- match.arg(dist)

  # 1) filter to one cut & threshold
  df <- subset(fuzzy, cut == cut & phi >= min_phi)
  if (nrow(df) == 0L) stop("No memberships at this cut/threshold.")

  # 2) build group x species matrix of phi (weights)
  #    xtabs keeps base-R only; rows=group, cols=species
  X <- xtabs(phi ~ group + species, data = df)
  X <- as.matrix(X)

  # Drop groups with all zeros (just in case)
  keep <- rowSums(X) > 0
  X <- X[keep, , drop = FALSE]
  if (nrow(X) < 2L) stop("Need >= 2 groups to cluster.")

  # 3) distance between groups
  if (dist == "bray") {
    # Bray-Curtis on nonnegative phi-weights (needs vegan)
    if (!requireNamespace("vegan", quietly = TRUE))
      stop("vegan package required for dist='bray'.")
    D <- vegan::vegdist(X, method = "bray")
  } else {
    # cosine distance (1 - cosine similarity), base-R
    row_norm <- sqrt(rowSums(X * X))
    row_norm[row_norm == 0] <- 1
    Xn <- X / row_norm
    S <- Xn %*% t(Xn)                    # cosine similarity
    S[S > 1] <- 1; S[S < -1] <- -1
    D <- 1 - S
    D <- stats::as.dist(D)
  }

  # 4) hclust
  hc <- stats::hclust(D, method = hclust_method)

  # 5) labels
  labs <- rownames(X)
  if (!is.null(group_names)) {
    # group_names should be a named vector: names = group codes
    labs <- ifelse(labs %in% names(group_names), group_names[labs], labs)
  }

  # 6) plot (optionally to PDF)
  if (!is.null(pdf_file)) {
    if (is.null(width_in)) {
      # ~0.25 inch per label, clamp 6-30"
      width_in <- min(30, max(6, 0.25 * length(labs)))
    }
    grDevices::pdf(pdf_file, width = width_in, height = height_in, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  plot(hc, labels = labs, main = sprintf("Fuzzy groups dendrogram (phi cut = %.2f, %s)", cut, dist),
       xlab = "", sub = "", cex = 0.8)

  invisible(list(hclust = hc, dist = D, X = X, labels = labs))
}
