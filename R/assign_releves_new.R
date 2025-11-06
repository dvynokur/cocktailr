#' Assign relevés (plots) to k groups from a releve_clustering result
#' with optional one-step child level (next branch division).
#'
#' @param x Object returned by `releve_clustering()` (class "releve_hclust").
#' @param k Integer, number of fine groups to cut the tree into (parent level).
#' @param next_level Logical; if TRUE, also compute the immediate children
#'   of each parent cluster (split each parent once into its two branches).
#'
#' @return A list with:
#'   - `groups`: data.frame(plot, cluster, cluster_orig, parent_cluster, child_cluster)
#'               (parent_cluster/child_cluster are included only if next_level=TRUE)
#'   - `sizes` : table of fine cluster sizes
#'   - `k`     : number of groups used
#'   - `hclust`: the hclust object (for convenience)
#'   - `child_map` : (when next_level=TRUE) integer vector mapping child_id -> parent_id
#' @export

assign_releves_new <- function(x, k, next_level = FALSE) {
  if (!inherits(x, "releve_hclust") || !"hclust" %in% names(x)) {
    stop("`x` must be the object returned by releve_clustering() (class 'releve_hclust').")
  }
  if (missing(k) || length(k) != 1L || !is.finite(k) || k < 1L) {
    stop("`k` must be a single positive integer.")
  }

  hc <- x$hclust
  # Use plot labels carried by releve_clustering if available; otherwise hclust labels
  labels <- if (!is.null(x$plot_labels)) x$plot_labels else hc$labels
  if (is.null(labels)) labels <- as.character(seq_along(hc$order))

  ## --- fine clusters at k (parent level) ---
  grp_orig <- stats::cutree(hc, k = k)        # named by hc$labels
  # Reorder to our `labels` vector for a stable row order
  if (!is.null(names(grp_orig))) grp_orig <- grp_orig[labels] else names(grp_orig) <- labels

  # Relabel clusters to 1..k by first leaf position along the dendrogram
  cl_ids <- sort(unique(grp_orig))
  ord <- hc$order                              # permutation of 1..n
  # positions of first leaf of each cluster along the dendrogram
  pos_first <- vapply(cl_ids, function(cl) {
    which_min <- which(ord %in% which(grp_orig == cl))[1]
    if (length(which_min) == 0L || is.na(which_min)) Inf else which_min
  }, integer(1))
  cl_order <- cl_ids[order(pos_first)]
  id_map <- stats::setNames(seq_along(cl_order), cl_order)
  grp_new <- unname(id_map[as.character(grp_orig)])

  # Base output (parent-only)
  out_df <- data.frame(
    plot = labels,
    cluster = grp_new,
    cluster_orig = grp_orig,
    stringsAsFactors = FALSE
  )

  sizes <- table(out_df$cluster)
  sizes <- sizes[order(as.integer(names(sizes)))]

  # If no child level requested, return here
  if (!isTRUE(next_level)) {
    return(list(
      groups = out_df,
      sizes  = sizes,
      k      = k,
      hclust = hc
    ))
  }

  ## --- next level: split each parent cluster once into its two branches ---
  # Convert to dendrogram for subtree operations
  dnd <- stats::as.dendrogram(hc)

  # Map label -> leaf order position (left-to-right)
  # hc$labels[hc$order] are the labels in plotting order
  ordered_labels <- hc$labels[hc$order]
  ord_pos <- stats::setNames(seq_along(ordered_labels), ordered_labels)

  # Helper: get all labels (leaves) underneath a dendrogram node
  get_leaves <- function(node) {
    if (is.leaf(node)) return(attr(node, "label"))
    c(get_leaves(node[[1]]), get_leaves(node[[2]]))
  }

  # Helper: find the subtree whose leaves exactly match `target_labels` (as a set)
  # Returns the node or NULL if not found.
  find_subtree_by_labels <- function(node, target_set) {
    labs <- get_leaves(node)
    if (length(labs) == length(target_set) && all(sort(labs) == sort(target_set))) {
      return(node)
    }
    if (is.leaf(node)) return(NULL)
    left  <- find_subtree_by_labels(node[[1]], target_set)
    if (!is.null(left))  return(left)
    right <- find_subtree_by_labels(node[[2]], target_set)
    if (!is.null(right)) return(right)
    NULL
  }

  parent_ids <- sort(unique(grp_new))   # 1..k in left-to-right order
  child_cluster <- integer(length(labels))  # will fill with child IDs
  child_map <- integer(0)                   # child_id -> parent_id
  next_child_id <- 0L

  for (pid in parent_ids) {
    # Labels (plots) belonging to this parent cluster
    lab_pid <- labels[grp_new == pid]

    if (length(lab_pid) <= 1L) {
      # Singleton parent: it can't split; make one child equal to the parent
      next_child_id <- next_child_id + 1L
      child_cluster[match(lab_pid, labels)] <- next_child_id
      child_map[next_child_id] <- pid
      next
    }

    # Find the subtree corresponding exactly to this parent's set of leaves
    subtree <- find_subtree_by_labels(dnd, lab_pid)
    if (is.null(subtree)) {
      # Fallback: if for any reason not found, assign as one child
      next_child_id <- next_child_id + 1L
      child_cluster[match(lab_pid, labels)] <- next_child_id
      child_map[next_child_id] <- pid
      next
    }

    # If the subtree is a leaf, same as singleton case
    if (is.leaf(subtree)) {
      next_child_id <- next_child_id + 1L
      child_cluster[match(lab_pid, labels)] <- next_child_id
      child_map[next_child_id] <- pid
      next
    }

    # Its two immediate children are the "next branch division"
    left_leaves  <- get_leaves(subtree[[1]])
    right_leaves <- get_leaves(subtree[[2]])

    # Order the two child branches left-to-right using first leaf position
    left_pos  <- min(ord_pos[left_leaves])
    right_pos <- min(ord_pos[right_leaves])

    if (left_pos <= right_pos) {
      first  <- left_leaves;  second <- right_leaves
    } else {
      first  <- right_leaves; second <- left_leaves
    }

    # Assign two child IDs
    next_child_id <- next_child_id + 1L
    child_cluster[match(first,  labels)] <- next_child_id
    child_map[next_child_id] <- pid

    next_child_id <- next_child_id + 1L
    child_cluster[match(second, labels)] <- next_child_id
    child_map[next_child_id] <- pid
  }

  # Add parent & child to output dataframe
  out_df$parent_cluster <- grp_new
  out_df$child_cluster  <- child_cluster

  # Recompute sizes for children as well (optional)
  # child_sizes <- table(out_df$child_cluster)

  list(
    groups    = out_df,
    sizes     = sizes,
    k         = k,
    hclust    = hc,
    child_map = child_map
  )
}
