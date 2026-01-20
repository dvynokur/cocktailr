#' List plots (relevés) belonging to Cocktail clusters or cluster unions
#'
#' @description
#' Returns the plot IDs (relevés) that belong to one or more Cocktail clusters
#' based on \code{x$Plot.cluster}. Cluster membership is defined as
#' \code{x$Plot.cluster > 0}.
#'
#' The \code{clusters} argument can define either:
#' \itemize{
#'   \item a vector of cluster IDs/labels (each element treated as its own cluster), or
#'   \item a list of vectors (each list element treated as a \strong{union group} of
#'         clusters; membership is the union across nodes).
#' }
#'
#' In union groups, a plot is considered a member if it is a member of
#' \emph{any} node in that group (OR logic).
#'
#' @param x A \code{"cocktail"} object (result of \code{\link{cocktail_cluster}}),
#'   containing \code{Plot.cluster}.
#'
#' @param clusters Cluster selection. Can be:
#' \itemize{
#'   \item a character vector of node labels like \code{c("c_12","c_27")}, or
#'   \item a numeric/integer vector like \code{c(12,27)}, or
#'   \item a list, where each element is a vector of labels/IDs defining a union group.
#' }
#'
#' @param values Logical; if \code{FALSE} (default), return plot IDs only.
#'   If \code{TRUE}, return per-plot membership \emph{values} for each cluster/group:
#'   \itemize{
#'     \item for a single cluster: the \code{Plot.cluster} values for that node;
#'     \item for a union group: the \emph{maximum} value across nodes (per plot),
#'           which corresponds to binary membership if \code{Plot.cluster} is binary,
#'           or the maximum relative cover value if \code{Plot.cluster} stores rel_cover.
#'   }
#'
#' @param drop0 Logical; if \code{TRUE} (default), drop plots with membership 0.
#'   If \code{FALSE}, keep all plots.
#'
#' @param return Character; return format:
#' \itemize{
#'   \item \code{"list"} (default): a list, one element per cluster/group.
#'   \item \code{"vector"}: only allowed when \code{clusters} defines exactly one
#'         cluster/group; returns a character vector of plot IDs (or numeric values if
#'         \code{values = TRUE}).
#'   \item \code{"matrix"}: return a matrix with plots in rows and clusters/groups
#'         in columns (only meaningful when \code{values = TRUE} or when keeping zeros).
#' }
#'
#' @return Depending on \code{return}:
#' \itemize{
#'   \item \code{"list"}: list of plot ID vectors (or named numeric vectors when
#'         \code{values = TRUE}).
#'   \item \code{"vector"}: plot IDs (or numeric membership values) for one group.
#'   \item \code{"matrix"}: numeric matrix with plots in rows and groups in columns.
#' }
#'
#' @export
releves_in_clusters <- function(
    x,
    clusters,
    values = FALSE,
    drop0  = TRUE,
    return = c("list", "vector", "matrix")
) {
  return <- match.arg(return)

  ## ---- checks -------------------------------------------------------------
  if (!is.list(x) || !"Plot.cluster" %in% names(x) || is.null(x$Plot.cluster)) {
    stop("`x` must be a Cocktail object with a non-NULL `Plot.cluster`.")
  }

  PC <- x$Plot.cluster
  if (!inherits(PC, "Matrix")) {
    PC <- Matrix::Matrix(as.matrix(PC), sparse = TRUE)
  }

  # Ensure column names exist and are in "c_<id>" form
  if (is.null(colnames(PC))) {
    colnames(PC) <- paste0("c_", seq_len(ncol(PC)))
  } else {
    num_only <- grepl("^\\d+$", colnames(PC))
    colnames(PC)[num_only] <- paste0("c_", colnames(PC)[num_only])
  }

  plot_names <- rownames(PC)
  if (is.null(plot_names)) {
    plot_names <- paste0("plot_", seq_len(nrow(PC)))
    rownames(PC) <- plot_names
  }

  if (missing(clusters) || is.null(clusters)) {
    stop("`clusters` must be provided (vector or list).")
  }

  ## ---- parse clusters into list-of-groups --------------------------------
  .parse_one <- function(v) {
    if (is.character(v)) {
      as.integer(sub("^c_", "", v))
    } else {
      as.integer(v)
    }
  }

  group_list <- if (is.list(clusters)) {
    lapply(clusters, .parse_one)
  } else {
    # vector -> each element is its own group
    lapply(as.list(clusters), .parse_one)
  }

  # clean groups
  n_nodes <- ncol(PC)
  group_list <- lapply(group_list, function(ids) {
    ids <- ids[is.finite(ids)]
    ids <- ids[ids > 0L & ids <= n_nodes]
    sort(unique(ids))
  })

  keep <- vapply(group_list, length, integer(1)) > 0L
  if (!any(keep)) stop("No valid cluster IDs found in `clusters`.")
  group_list <- group_list[keep]

  ## ---- compute membership -------------------------------------------------
  out_vals <- vector("list", length(group_list))
  out_lab  <- character(length(group_list))

  for (g in seq_along(group_list)) {
    ids  <- group_list[[g]]
    cols <- paste0("c_", ids)

    missing_cols <- setdiff(cols, colnames(PC))
    if (length(missing_cols)) {
      stop("`x$Plot.cluster` is missing columns: ", paste(missing_cols, collapse = ", "))
    }

    sub_pc <- PC[, cols, drop = FALSE]

    if (values) {
      # union group: max value per plot
      v <- if (ncol(sub_pc) == 1L) {
        as.numeric(sub_pc[, 1])
      } else {
        apply(as.matrix(sub_pc), 1, max)
      }
      names(v) <- plot_names
      if (drop0) v <- v[v > 0]
      out_vals[[g]] <- v

    } else {
      # membership only (binary)
      m <- Matrix::rowSums(sub_pc != 0) > 0
      plots <- plot_names[m]
      out_vals[[g]] <- plots
    }

    # group label
    out_lab[g] <- if (length(ids) == 1L) {
      paste0("c_", ids)
    } else {
      paste0("g_", paste(ids, collapse = "_"))
    }
  }

  names(out_vals) <- out_lab

  ## ---- return formats -----------------------------------------------------
  if (return == "vector") {
    if (length(out_vals) != 1L) {
      stop("`return = \"vector\"` requires exactly one cluster/group.")
    }
    return(out_vals[[1L]])
  }

  if (return == "matrix") {
    # meaningful mainly when values=TRUE or drop0=FALSE
    M <- matrix(0, nrow = length(plot_names), ncol = length(out_vals),
                dimnames = list(plot_names, names(out_vals)))

    for (j in seq_along(out_vals)) {
      if (values) {
        v <- out_vals[[j]]
        M[names(v), j] <- v
      } else {
        M[out_vals[[j]], j] <- 1
      }
    }

    if (drop0) {
      keep_rows <- rowSums(M != 0) > 0
      M <- M[keep_rows, , drop = FALSE]
    }
    return(M)
  }

  # return == "list"
  out_vals
}
