#' Save cocktail clustering result
#'
#' Export the result of cocktail_cluster() to CSV, RDS, or both.
#'
#' @param res A cocktail result object from cocktail_cluster().
#' @param outdir Directory for CSVs (created if missing) when \code{format} is
#'        "csv" or "both".
#' @param format One of \code{"csv"}, \code{"rds"}, \code{"both"}.
#' @param rds_file Optional file path for the RDS. If NULL and \code{outdir} is
#'        provided, defaults to \file{<outdir>/cocktail_result.rds}.
#' @param rds_compress Compression for RDS: \code{"xz"}, \code{"gzip"},
#'        \code{"bzip2"}, or \code{"none"}.
#'
#' @return Invisibly TRUE on success.
#' @importFrom utils write.csv
#' @export
cocktail_save <- function(res,
                          outdir = NULL,
                          format = c("csv","rds","both"),
                          rds_file = NULL,
                          rds_compress = c("xz","gzip","bzip2","none")) {

  format       <- match.arg(format)
  rds_compress <- match.arg(rds_compress)

  if (format %in% c("csv","both")) {
    if (is.null(outdir)) stop("Please provide `outdir` for CSV output.")
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    utils::write.csv(res$Cluster.species, file.path(outdir,"Cluster_species.csv"), row.names = FALSE)
    utils::write.csv(res$Cluster.info,    file.path(outdir,"Cluster_info.csv"),    row.names = FALSE)
    utils::write.csv(res$Plot.cluster,    file.path(outdir,"Plot_cluster.csv"),    row.names = FALSE)
    utils::write.csv(res$Cluster.merged,  file.path(outdir,"Cluster_merged.csv"),  row.names = FALSE)
    utils::write.csv(data.frame(height = res$Cluster.height),
                     file.path(outdir,"Cluster_height.csv"), row.names = FALSE)
    utils::write.csv(data.frame(species = res$species),
                     file.path(outdir,"species.csv"), row.names = FALSE)
    utils::write.csv(data.frame(plots   = res$plots),
                     file.path(outdir,"plots.csv"),   row.names = FALSE)
  }

  if (format %in% c("rds","both")) {
    if (is.null(rds_file)) {
      if (is.null(outdir)) stop("Provide `rds_file` or `outdir` for RDS.")
      rds_file <- file.path(outdir, "cocktail_result.rds")
    }
    saveRDS(res, rds_file, compress = rds_compress)
  }

  invisible(TRUE)
}
