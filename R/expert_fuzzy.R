#' Export fuzzy Cocktail groups to JUICE expert system format
#'
#' @param fuzzy_result Data.frame from cocktail_fuzzy(),
#'   with columns: species, group, phi.
#' @param phi_min Minimum fidelity (phi) for a species to be included
#'   in a group (default 0.3).
#' @param weighted Logical. If FALSE (default), species are listed
#'   alphabetically. If TRUE, species are sorted by decreasing phi
#'   and phi*100 is prefixed.
#' @param file Output path for the expert system text file.
#'
#' @return Invisibly returns the filtered table used for export.
#' @export
juice_expert_system <- function(fuzzy_result,
                                phi_min = 0.3,
                                weighted = FALSE,
                                file = "expert_system.txt") {
  stopifnot(all(c("species", "group", "phi") %in% names(fuzzy_result)))

  # Filter species by phi threshold
  dat <- subset(fuzzy_result, phi >= phi_min)
  if (nrow(dat) == 0) {
    stop("No species above phi threshold = ", phi_min)
  }

  # Split by group
  groups <- split(dat, dat$group)

  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  for (g in names(groups)) {
    gdat <- groups[[g]]

    # Group name = group code + " - " + top species
    top_spp <- gdat$species[which.max(gdat$phi)]
    header <- paste0("### ", g, " (", top_spp, " group)")
    writeLines(header, con)

    if (weighted) {
      # Sort by decreasing phi
      gdat <- gdat[order(-gdat$phi, gdat$species), ]
      lines <- sprintf("%5.1f %s", 100 * gdat$phi, gdat$species)
    } else {
      # Sort alphabetically
      gdat <- gdat[order(gdat$species), ]
      lines <- paste0("     ", gdat$species)
    }
    writeLines(lines, con)
    writeLines("", con) # blank line after each block
  }

  invisible(dat)
}
