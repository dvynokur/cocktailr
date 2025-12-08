cocktailr
================

- [cocktailr](#cocktailr)
  - [Overview](#overview)
  - [Background](#background)
  - [Installation](#installation)
  - [Quick Start](#quick-start)
  - [Typical workflow](#typical-workflow)
    - [1. Visualise the dendrogram](#1-visualise-the-dendrogram)
    - [2. Select parent clusters at a φ
      cut](#2-select-parent-clusters-at-a-φ-cut)
    - [3. Diagnostic species for parent
      clusters](#3-diagnostic-species-for-parent-clusters)
    - [4. φ-based distances between
      clusters](#4-φ-based-distances-between-clusters)
    - [5. Assign plots (relevés) to
      groups](#5-assign-plots-relevés-to-groups)
    - [(Optional) Attach assignments to a header data
      frame](#optional-attach-assignments-to-a-header-data-frame)
  - [Reference](#reference)

# cocktailr

Fast, reproducible *Cocktail* clustering for vegetation tables.

------------------------------------------------------------------------

## Overview

**cocktailr** provides fast and reproducible *Cocktail* clustering of
vegetation data, identifying groups of co-occurring species from **plots
× species** tables. It uses optimized sparse-matrix calculations and φ
(phi) coefficients to produce consistent, deterministic results, even
for large vegetation databases.

The package implements:

- **Hierarchical Cocktail clustering** of species
  (`cocktail_cluster()`).
- **Dendrogram plotting** with φ heights and optional cluster bands
  (`cocktail_plot()`).
- Extraction of **parent clusters at a φ cut** (`clusters_at_cut()`).
- **Diagnostic species lists** for clusters or unions of clusters
  (`species_in_clusters()`).
- **φ-based distances between clusters** (`cluster_phi_dist()`).
- **Assignment of plots (relevés) to groups** using several strategies
  that combine covers, topology, and species–cluster φ
  (`assign_releves()`).

Fuzzy / graded affinities are handled via the optional **species ×
cluster φ matrix** (`Species.cluster.phi`) that `cocktail_cluster()` can
compute when `species_cluster_phi = TRUE`.

------------------------------------------------------------------------

## Background

The *Cocktail* method (Bruelheide 2000, 2016) identifies sets of species
that co-occur more often than expected by chance and merges them
hierarchically according to the **phi coefficient of association**. Each
resulting cluster is characterized by its diagnostic species and a
threshold (*m*) indicating how many group species a plot must contain to
belong to it.

For details, see the original works:

- Bruelheide, H. (2000). *A new measure of fidelity and its application
  to defining species groups.* **Journal of Vegetation Science**, 11,
  167–178. <https://doi.org/10.2307/3236796>  
- Bruelheide, H. (2016). *Cocktail clustering – a new hierarchical
  agglomerative algorithm for extracting species groups in vegetation
  databases.* **Journal of Vegetation Science**, 27(6), 1297–1307.
  <https://doi.org/10.1111/jvs.12454>

------------------------------------------------------------------------

## Installation

``` r
# Install from GitHub
remotes::install_github("dvynokur/cocktailr")
```

------------------------------------------------------------------------

## Quick Start

A minimal example that demonstrates a merge with **positive φ** and
shows non-empty cluster results.

``` r
library(cocktailr)

# Tiny matrix with positive association between sp1 & sp2
vm <- matrix(
  c(1, 1,
    1, 0,
    0, 0),
  nrow = 3, byrow = TRUE,
  dimnames = list(paste0("plot", 1:3), c("sp1", "sp2"))
)

# Basic Cocktail clustering
res <- cocktail_cluster(vm, progress = FALSE)

# Merge height (phi) for the only internal node
res$Cluster.height
#> [1] 0.5

# Parent clusters at phi >= 0.3
labs <- clusters_at_cut(res, phi = 0.3)
labs
#> [1] "c_1"

# Topological species sets for those clusters
species_in_clusters(res, labels = labs)
#> $c_1
#> [1] "sp1" "sp2"
```

------------------------------------------------------------------------

## Typical workflow

A small end-to-end example on a toy **plots × species** matrix, showing:

1.  Cocktail clustering with optional species–cluster φ.
2.  Dendrogram plotting.
3.  Selecting clusters at a φ cut.
4.  Diagnostic species lists.
5.  φ-based distances between clusters and grouping.
6.  Plot assignment using φ and cover.

``` r
library(cocktailr)

# Toy plots × species matrix
vm <- matrix(
  c(1, 0, 1, 0, 1, 0,
    0, 1, 0, 1, 0, 1,
    1, 1, 0, 0, 1, 0,
    0, 0, 1, 1, 0, 0),
  nrow = 4, byrow = TRUE,
  dimnames = list(
    paste0("plot", 1:4),
    paste0("sp",   1:6)
  )
)

# 1) Cocktail clustering, keeping relative cover and species–cluster phi
res <- cocktail_cluster(
  vegmatrix           = vm,
  progress            = FALSE,
  plot_values         = "rel_cover",     # keeps cover-based Plot.cluster
  species_cluster_phi = TRUE             # computes Species.cluster.phi
)

names(res)
#> [1] "Cluster.species"     "Cluster.info"        "Plot.cluster"       
#> [4] "Cluster.merged"      "Cluster.height"      "Species.cluster.phi"
#> [7] "species"             "plots"
```

### 1. Visualise the dendrogram

``` r
# Plot to the current device with a phi cut and labels at that cut
cocktail_plot(
  x              = res,
  file           = NULL,       # RStudio Plots pane / current device
  phi_cut        = 0.3,
  label_clusters = TRUE,
  cex_species    = 0.9
)
```

<img src="man/figures/README-typical-plot-1.png" width="100%" />

### 2. Select parent clusters at a φ cut

``` r
phi_cut <- 0.3

parent_labels <- clusters_at_cut(
  x         = res,
  phi       = phi_cut,
  as_labels = TRUE
)

parent_labels
#> [1] "c_1" "c_3"
```

### 3. Diagnostic species for parent clusters

Topological species per parent cluster:

``` r
diag_sp_topo <- species_in_clusters(
  x      = res,
  labels = parent_labels
)

diag_sp_topo
#> $c_1
#> [1] "sp1" "sp5"
#> 
#> $c_3
#> [1] "sp2" "sp4" "sp6"
```

With φ-based filtering and ranking (uses `Species.cluster.phi`):

``` r
diag_sp_phi <- species_in_clusters(
  x                   = res,
  labels              = parent_labels,
  species_cluster_phi = TRUE,
  min_phi             = 0.20,
  top_k               = 10
)

diag_sp_phi
#> $c_1
#>   species phi
#> 1     sp1   1
#> 2     sp5   1
#> 
#> $c_3
#>   species       phi
#> 1     sp6 1.0000000
#> 2     sp2 0.5773503
#> 3     sp4 0.5773503
```

### 4. φ-based distances between clusters

Compute a distance matrix between clusters based on species fidelity
profiles:

``` r
D <- cluster_phi_dist(
  x        = res,
  clusters = parent_labels   # could also be node IDs, e.g. c(5, 12, 18)
)

D

# Example: hierarchical clustering of clusters
hc_nodes <- hclust(D, method = "average")
plot(hc_nodes, main = "Cluster dendrogram (phi-based distance)")
```

You can inspect a similar tree in your own data and decide to merge
similar parent clusters if desired.

### 5. Assign plots (relevés) to groups

Use `assign_releves()` with one of the strategies:

- `"count"` – number of diagnostic species present.
- `"cover"` – summed cover of diagnostic species.
- `"phi_topo"` – sum of phi over topological species present.
- `"phi_cover_topo"` – sum of cover × phi over topological species.
- `"phi_cover"` – sum of cover × phi for species with phi ≥ `min_phi`.
- `"phi"` – sum of phi for species with phi ≥ `min_phi`.

Here we use `"phi_cover"` with the parent clusters at `phi_cut`:

``` r
assign_phi <- assign_releves(
  x              = res,
  vegmatrix      = vm,
  strategy       = "phi_cover",
  phi_cut        = phi_cut,
  min_phi        = 0.20,
  min_group_size = 1L
)

assign_phi
#> plot1 plot2 plot3 plot4 
#> "g_1" "g_3" "g_1"    NA 
#> attr(,"details")
#> attr(,"details")$strategy
#> [1] "phi_cover"
#> 
#> attr(,"details")$phi_cut
#> [1] 0.3
#> 
#> attr(,"details")$clusters
#> NULL
#> 
#> attr(,"details")$groups_used
#> [1] "g_1" "g_3"
#> 
#> attr(,"details")$min_phi
#> [1] 0.2
#> 
#> attr(,"details")$min_group_size
#> [1] 1
#> 
#> attr(,"details")$collapsed_groups
#> character(0)
table(assign_phi)
#> assign_phi
#> g_1 g_3 
#>   2   1
```

The returned vector is named by plot ID and contains:

- group labels like `"g_5"` (or combinations if you defined union
  groups);
- `"+"` for ties between groups;
- `"-"` for groups that were collapsed by `min_group_size`;
- `NA` when no group wins according to the chosen strategy.

You can then attach these assignments to your plot header / metadata
table.

------------------------------------------------------------------------

### (Optional) Attach assignments to a header data frame

If you have a header table `hea` with a column `releve_number`, you can
add multiple assignment strategies as new columns:

``` r
library(dplyr)

# example strategies you want as separate columns
strategies <- c("count", "cover", "phi_topo", "phi_cover_topo", "phi_cover")

hea2 <- hea

for (s in strategies) {
  rel_assigned <- assign_releves(
    x              = res,
    vegmatrix      = vm,
    strategy       = s,
    phi_cut        = 0.3,
    min_phi        = 0.2,
    min_group_size = 2
  )

  # align by releve_number using the names of rel_assigned
  idx <- match(hea2$releve_number, as.integer(names(rel_assigned)))

  colname <- paste0("grp_", s)  # e.g. "grp_count", "grp_cover", ...

  hea2[[colname]] <- rel_assigned[idx]

  # optional: turn "-" into NA instead of a literal dash
  # hea2[[colname]][hea2[[colname]] == "-"] <- NA_character_
}

dplyr::glimpse(hea2)
```

------------------------------------------------------------------------

## Reference

See function help for details:

- `?cocktail_cluster` – build the Cocktail tree, optionally with
  species–cluster phi  
- `?cocktail_plot` – draw dendrograms (PDF/PNG or current device)  
- `?clusters_at_cut` – parent clusters at a phi cut  
- `?species_in_clusters` – diagnostic species per node or node union  
- `?cluster_phi_dist` – phi-based distances between clusters  
- `?assign_releves` – assign plots to groups using covers and phi

------------------------------------------------------------------------

© 2025 Denys Vynokurov & Helge Bruelheide. Licensed under MIT.
