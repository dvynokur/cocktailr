#' Find Cocktail clusters that contain given species
#'
#' @description
#' Returns the internal Cocktail nodes (clusters) whose node-constituting species
#' set (from \code{x$Cluster.species}) includes the queried species.
#'
#' This function uses \strong{topological membership} only (species in node as stored
#' in \code{x$Cluster.species}), not fidelity weights. It is therefore stable and fast.
#'
#' By default, the function excludes clusters with merge height \eqn{\phi < 0}
#' (via \code{min_phi = 0}), which typically removes very weak / root-like merges.
#'
#' @param x A \code{"cocktail"} object (result of \code{\link{cocktail_cluster}}),
#'   containing \code{Cluster.species}. Components \code{Cluster.height} and
#'   \code{Cluster.info} (with columns \code{"k","m"}) are used for filtering and
#'   table output.
#'
#' @param species Character vector of species names to query. Must match column
#'   names of \code{x$Cluster.species}. Species not present in the Cocktail object
#'   are ignored with a warning (unless none match).
#'
#' @param match Character; how to match multiple species:
#' \itemize{
#'   \item \code{"any"} (default): return clusters containing \emph{at least one}
#'         of the requested species.
#'   \item \code{"all"}: return clusters containing \emph{all} requested species.
#' }
#'
#' @param min_phi Numeric scalar in \eqn{[-1,1]} (default \code{0}).
#'   Only clusters with \code{x$Cluster.height >= min_phi} are returned.
#'
#' @param return Character; what to return:
#' \itemize{
#'   \item \code{"labels"} (default): character labels like \code{"c_12"}.
#'   \item \code{"ids"}: integer node IDs.
#'   \item \code{"table"}: data frame of selected nodes, including
#'     \code{h,k,m} if available.
#' }
#'
#' @return Depending on \code{return}:
#' \itemize{
#'   \item \code{"labels"}: character vector of node labels.
#'   \item \code{"ids"}: integer vector of node IDs.
#'   \item \code{"table"}: data frame with columns \code{cluster, n_match, h, k, m}.
#' }
#'
#' @export
clusters_with_species <- function(
    x,
    species,
    match   = c("any", "all"),
    min_phi = 0,
    return  = c("labels", "ids", "table")
) {
  match  <- match.arg(match)
  return <- match.arg(return)

  ## ---- checks -------------------------------------------------------------
  if (!is.list(x) || !"Cluster.species" %in% names(x) || is.null(x$Cluster.species)) {
    stop("`x` must be a Cocktail object with a non-NULL `Cluster.species`.")
  }

  CS <- x$Cluster.species
  if (!is.matrix(CS)) stop("`x$Cluster.species` must be a matrix.")
  if (is.null(colnames(CS))) stop("`x$Cluster.species` must have species column names.")

  if (missing(species) || is.null(species) || !length(species)) {
    stop("`species` must be a non-empty character vector.")
  }
  species <- as.character(species)

  # min_phi validity
  if (!is.numeric(min_phi) || length(min_phi) != 1L || is.na(min_phi)) {
    stop("`min_phi` must be a single numeric value in [-1,1].")
  }
  if (min_phi < -1 || min_phi > 1) {
    stop("`min_phi` must be a numeric value in [-1,1].")
  }
  if (is.null(x$Cluster.height) || length(x$Cluster.height) < nrow(CS)) {
    stop("Filtering by `min_phi` requires `x$Cluster.height` (from cocktail_cluster()).")
  }

  ## ---- validate species ---------------------------------------------------
  spp_all <- colnames(CS)
  present <- intersect(species, spp_all)
  missing <- setdiff(species, spp_all)

  if (length(missing)) {
    warning(
      "Some species are not present in `x$Cluster.species` and will be ignored: ",
      paste(head(missing, 10), collapse = ", "),
      if (length(missing) > 10) " ..." else ""
    )
  }

  if (!length(present)) {
    stop("None of the requested species are present in `x$Cluster.species`.")
  }

  ## ---- compute cluster hits ----------------------------------------------
  sub <- CS[, present, drop = FALSE] > 0

  if (match == "any") {
    hit <- rowSums(sub) > 0
  } else {
    hit <- rowSums(sub) == length(present)
  }

  ids <- which(hit)
  if (!length(ids)) {
    if (return == "table") {
      return(data.frame(
        cluster = character(0),
        n_match = integer(0),
        h       = numeric(0),
        k       = numeric(0),
        m       = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    if (return == "labels") return(character(0))
    return(integer(0))
  }

  ## ---- phi filtering (default min_phi = 0) -------------------------------
  H <- x$Cluster.height
  ids <- ids[H[ids] >= min_phi]

  if (!length(ids)) {
    if (return == "table") {
      return(data.frame(
        cluster = character(0),
        n_match = integer(0),
        h       = numeric(0),
        k       = numeric(0),
        m       = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    if (return == "labels") return(character(0))
    return(integer(0))
  }

  ## ---- return ids/labels --------------------------------------------------
  if (return == "ids") {
    return(ids)
  }

  if (return == "labels") {
    return(paste0("c_", ids))
  }

  ## ---- table output -------------------------------------------------------
  n_match <- rowSums(sub[ids, , drop = FALSE])

  # optional fields
  h <- as.numeric(H[ids])

  if (!is.null(x$Cluster.info) && all(c("k", "m") %in% colnames(x$Cluster.info))) {
    k <- as.numeric(x$Cluster.info[ids, "k"])
    m <- as.numeric(x$Cluster.info[ids, "m"])
  } else {
    k <- rep(NA_real_, length(ids))
    m <- rep(NA_real_, length(ids))
  }

  data.frame(
    cluster = paste0("c_", ids),
    n_match = as.integer(n_match),
    h       = h,
    k       = k,
    m       = m,
    stringsAsFactors = FALSE
  )
}
