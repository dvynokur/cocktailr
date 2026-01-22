#' cocktailr: fast Cocktail clustering for vegetation tables
#'
#' **cocktailr** provides fast and reproducible *Cocktail* clustering of vegetation data,
#' identifying groups of co-occurring species from **plots × species** tables.
#' It uses sparse-matrix calculations and φ (phi) coefficients to produce deterministic
#' results, even for large vegetation databases.
#'
#' @section Overview:
#' The package implements:
#' \itemize{
#'   \item \code{\link{cocktail_cluster}} — hierarchical Cocktail clustering of species.
#'   \item \code{\link{cocktail_plot}} — dendrogram plotting with φ heights and optional cluster bands.
#'   \item \code{\link{clusters_at_cut}} — extract parent clusters at a φ cut.
#'   \item \code{\link{select_clusters}} — select strong clusters by a combined score.
#'   \item \code{\link{species_in_clusters}} — diagnostic species lists for clusters or cluster unions.
#'   \item \code{\link{releves_in_clusters}} — list plots (relevés) with membership in clusters or cluster unions.
#'   \item \code{\link{clusters_with_species}} — find clusters associated with a given species (or set of species).
#'   \item \code{\link{cluster_phi_dist}} — distances between clusters based on plot co-membership φ.
#'   \item \code{\link{assign_releves}} — assign plots (relevés) to groups using cover- and φ-based strategies.
#' }
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item Run clustering: \code{\link{cocktail_cluster}}
#'   \item Visualise dendrogram: \code{\link{cocktail_plot}}
#'   \item Select clusters: \code{\link{clusters_at_cut}} or \code{\link{select_clusters}}
#'   \item Diagnose clusters: \code{\link{species_in_clusters}}, \code{\link{releves_in_clusters}},
#'         \code{\link{clusters_with_species}}
#'   \item Compare clusters (optional): \code{\link{cluster_phi_dist}}
#'   \item Assign plots to groups: \code{\link{assign_releves}}
#' }
#'
#' @section Background:
#' The *Cocktail* method (Bruelheide 2000, 2016) merges species into groups based on the
#' **phi coefficient of association**. Each cluster is characterised by its member species
#' and a membership threshold (*m*) indicating how many cluster species a plot must contain
#' to belong to that node.
#'
#' @section References:
#' \itemize{
#'   \item Bruelheide, H. (2000). A new measure of fidelity and its application to defining species groups.
#'         \emph{Journal of Vegetation Science}, 11, 167–178. \doi{10.2307/3236796}
#'   \item Bruelheide, H. (2016). Cocktail clustering – a new hierarchical agglomerative algorithm for extracting
#'         species groups in vegetation databases. \emph{Journal of Vegetation Science}, 27(6), 1297–1307.
#'         \doi{10.1111/jvs.12454}
#' }
#'
#' @seealso
#' \itemize{
#'   \item GitHub: \url{https://github.com/dvynokur/cocktailr}
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import Matrix
## usethis namespace: end
NULL
