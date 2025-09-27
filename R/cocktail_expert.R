#' Build a JUICE-readable expert system from Cocktail results
#'
#' @description
#' Produce a plain-text expert system where each group starts with a header line
#' `"### <group-name>"` followed by one species per line.
#' Works with either a Cocktail object (species lists per cluster) or a fuzzy
#' \eqn{\phi} matrix (species \u00D7 nodes).
#'
#' @param x Either:
#'   \itemize{
#'     \item a Cocktail object (list with `Cluster.species`), or
#'     \item a numeric matrix returned by `cocktail_fuzzy()` (species \u00D7 nodes;
#'           columns named like `"c_2127"`).
#'   }
#' @param labels Character labels like `"c_2127"` or integer node IDs to include.
#'   If `x` is a \eqn{\phi} matrix and `labels` is `NULL`, all columns are used.
#'   If `x` is a Cocktail object, `labels` must be supplied.
#' @param min_phi Fuzzy only: keep species with \eqn{\phi \ge} `min_phi`. Default `0.20`.
#' @param top_k Fuzzy only: optional integer; after `min_phi` filtering keep only the
#'   top `k` species by \eqn{\phi}. Default `NULL` = keep all that pass `min_phi`.
#' @param fuzzy_sort Fuzzy only: how to order species in the output —
#'   `"phi"` (default; decreasing \eqn{\phi}), `"alpha"` (alphabetical by species),
#'   or `"none"` (input order). **Note:** group names are always derived from the
#'   highest-\eqn{\phi} species, regardless of this setting.
#' @param group_naming Fuzzy only: how to name groups in the `"###"` header —
#'   `"id"` (use the column label, e.g., `"c_2127"`), `"top1"` (highest-\eqn{\phi}
#'   species + `-group`), or `"top2"` (top two \eqn{\phi} species joined by `-`,
#'   plus `-group`). Spaces are replaced with `-` and non-alphanumeric characters
#'   are removed.
#' @param file Optional path to write a `.txt`. If `NULL` (default), the function
#'   returns the text as a single character string.
#'
#' @details
#' For Cocktail objects, each species is printed on its own line, indented with
#' five spaces. For \eqn{\phi} matrices, each line starts with `100 * phi`
#' formatted to one decimal place (e.g., `44.6`), followed by a space and the
#' species name. A blank line separates groups.
#'
#' @return A single character string containing the expert-system text. If
#'   `file` is provided, the text is written to disk and returned invisibly.
#'
#' @seealso `species_in_clusters()`, `cocktail_fuzzy()`, `clusters_at_cut()`
#' @export

cocktail_expert <- function(
    x,
    labels       = NULL,
    min_phi      = 0.20,
    top_k        = NULL,
    fuzzy_sort   = c("phi", "alpha", "none"),
    group_naming = c("id", "top1", "top2"),
    file         = NULL
) {
  fuzzy_sort   <- match.arg(fuzzy_sort)
  group_naming <- match.arg(group_naming)

  .is_cocktail <- function(obj) is.list(obj) && "Cluster.species" %in% names(obj)

  .norm_labels_to_chars <- function(lab, available) {
    if (is.null(lab)) return(NULL)
    if (is.character(lab)) {
      cand <- ifelse(grepl("^c_\\d+$", lab), lab, paste0("c_", sub("^c_", "", lab)))
    } else if (is.numeric(lab) || is.integer(lab)) {
      cand <- paste0("c_", as.integer(lab))
    } else stop("`labels` must be character like 'c_123' or integer node IDs.")
    hit <- intersect(cand, available)
    miss <- setdiff(cand, hit)
    if (length(miss)) warning("Dropping unknown labels: ", paste(miss, collapse = ", "))
    hit
  }

  .sanitize_group_name <- function(s) {
    s <- gsub("[[:space:]]+", "-", s)
    s <- gsub("[^-[:alnum:]_]", "", s)
    s
  }

  out_blocks <- list()

  if (.is_cocktail(x)) {
    ## ---- Cocktail object: plain species lists ----
    CS <- x$Cluster.species
    sp_names <- colnames(CS)
    if (is.null(sp_names)) stop("Cluster.species must have column (species) names.")
    if (is.null(labels)) stop("For Cocktail objects, please supply `labels`.")

    available <- paste0("c_", seq_len(nrow(CS)))
    labs_use  <- .norm_labels_to_chars(labels, available)
    if (!length(labs_use)) return("")

    for (lab in labs_use) {
      node_id <- as.integer(sub("^c_", "", lab))
      spp <- sp_names[which(CS[node_id, ] == 1L)]
      block <- c(
        paste0("### ", .sanitize_group_name(lab)),
        if (length(spp)) paste0("     ", spp) else character(0),
        ""  # blank line separator
      )
      out_blocks[[length(out_blocks) + 1L]] <- block
    }

  } else {
    ## ---- Fuzzy φ matrix: value + species per group ----
    if (!is.matrix(x) || !is.numeric(x)) {
      stop("`x` must be either a Cocktail object or a numeric φ matrix from cocktail_fuzzy().")
    }
    Phi <- x
    if (is.null(rownames(Phi))) stop("φ matrix must have species in rownames.")
    if (is.null(colnames(Phi))) stop("φ matrix must have node labels in colnames (e.g., 'c_2127').")

    labs_use <- if (is.null(labels)) colnames(Phi) else .norm_labels_to_chars(labels, colnames(Phi))
    if (!length(labs_use)) return("")

    if (!is.null(top_k)) {
      top_k <- as.integer(top_k)
      if (is.na(top_k) || top_k <= 0) top_k <- NULL
    }
    min_phi <- as.numeric(min_phi); if (!is.finite(min_phi)) min_phi <- 0

    for (lab in labs_use) {
      v  <- as.numeric(Phi[, lab])
      sp <- rownames(Phi)

      # ----- Group naming always based on top-φ species (before filtering or display sorting)
      ord_phi <- order(v, decreasing = TRUE, na.last = NA)
      top1 <- if (length(ord_phi) >= 1) sp[ord_phi[1]] else ""
      top2 <- if (length(ord_phi) >= 2) sp[ord_phi[2]] else ""
      gname <- lab
      if (group_naming == "top1" && nzchar(top1)) {
        gname <- paste0(.sanitize_group_name(top1), "-group")
      } else if (group_naming == "top2" && nzchar(top1)) {
        nm <- if (nzchar(top2)) paste(top1, top2, sep = "-") else top1
        gname <- paste0(.sanitize_group_name(nm), "-group")
      }
      gname <- .sanitize_group_name(gname)

      # ----- Display table (filter/sort as requested)
      df <- data.frame(species = sp, phi = v, stringsAsFactors = FALSE)
      df <- df[is.finite(df$phi), , drop = FALSE]
      df <- df[df$phi >= min_phi, , drop = FALSE]

      if (fuzzy_sort == "phi") {
        df <- df[order(df$phi, decreasing = TRUE), , drop = FALSE]
      } else if (fuzzy_sort == "alpha") {
        df <- df[order(df$species, decreasing = FALSE), , drop = FALSE]
      } # else "none"

      if (!is.null(top_k) && nrow(df) > top_k) df <- df[seq_len(top_k), , drop = FALSE]

      # lines like: "44.6 Species name"
      lines <- if (nrow(df)) sprintf("%.1f %s", 100 * df$phi, df$species) else character(0)

      block <- c(
        paste0("### ", gname),
        lines,
        ""  # blank line separator
      )
      out_blocks[[length(out_blocks) + 1L]] <- block
    }
  }

  txt <- paste(unlist(out_blocks, use.names = FALSE), collapse = "\n")
  if (!is.null(file)) {
    writeLines(txt, file, useBytes = TRUE)
    invisible(txt)
  } else {
    txt
  }
}
