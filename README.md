cocktailr
================

- [cocktailr](#cocktailr)
  - [Overview](#overview)
  - [Background](#background)
  - [Installation](#installation)
  - [Quick Start](#quick-start)
  - [Typical workflow](#typical-workflow)
    - [(Optional) Attach assignments to header data
      frame](#optional-attach-assignments-to-header-data-frame)
  - [Reference](#reference)

# cocktailr

Fast, reproducible *Cocktail* clustering for vegetation tables.

------------------------------------------------------------------------

## Overview

**cocktailr** provides fast and reproducible *Cocktail* clustering of
vegetation data, identifying groups of co-occurring species from **plots
× species** tables. It uses optimized sparse-matrix calculations and
exact φ (phi) coefficients to produce consistent, deterministic results.
The package can also estimate fuzzy memberships, showing how strongly
each species is associated with different clusters.

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
vm <- matrix(c(
  1,1,
  1,0,
  0,0
), nrow = 3, byrow = TRUE,
dimnames = list(paste0("plot", 1:3), c("sp1","sp2")))

res  <- cocktail_cluster(vm, progress = FALSE)
res$Cluster.height  # should be > 0 (e.g. +0.5)
#> [1] 0.5

# Parent clusters at φ ≥ 0.3
labs <- clusters_at_cut(res, phi = 0.3)
labs
#> [1] "c_1"
species_in_clusters(res, labels = labs)
#> $c_1
#> [1] "sp1" "sp2"
```

------------------------------------------------------------------------

## Typical workflow

An end-to-end example on a toy **plots × species** matrix.  
It shows classical clustering, dendrogram plotting, cluster inspection,
fuzzy φ computation, and plot assignment.

``` r
library(cocktailr)

# Toy plots × species matrix
vm <- matrix(c(1,0,1,0,1,0,
               0,1,0,1,0,1,
               1,1,0,0,1,0,
               0,0,1,1,0,0),
             nrow = 4, byrow = TRUE,
             dimnames = list(paste0("plot", 1:4), paste0("sp", 1:6)))

# 1) Classical Cocktail clustering
res <- cocktail_cluster(vm, progress = FALSE)

# 2) Plot dendrogram to PDF (disabled in README to avoid file output)
# plot_cocktail(res, "res.pdf", bands_phi = 0.30, palette = "rainbow", cex_labels = 1)

# 3) Parent clusters at φ ≥ 0.30
phi_cut <- 0.30
labs <- clusters_at_cut(res, phi = phi_cut)
labs
#> [1] "c_1" "c_3"

# 4) Diagnostic species from Cocktail object (binary membership)
species_in_clusters(res, labels = labs)
#> $c_1
#> [1] "sp1" "sp5"
#> 
#> $c_3
#> [1] "sp2" "sp4" "sp6"

# 5) Fuzzy Cocktail: φ for species × all internal nodes
res_fuzzy <- cocktail_fuzzy(res, vm)

# 6) Diagnostic species from fuzzy φ matrix
species_in_clusters(res_fuzzy, labels = labs, min_phi = 0.20, top_k = 10)
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

# 7) Assign plots to parent groups (strict)
assign_strict <- assign_releves(
  x               = res,
  vegmatrix       = vm,
  mode            = "strict",
  phi             = phi_cut,
  cover_transform = "sqrt",
  compare         = "max",
  min_cover       = 0,
  min_group_size  = 1
)
table(assign_strict)
#> assign_strict
#> c_1 c_3 
#>   2   2

# 8) Assign plots using fuzzy φ weights
assign_fuzzy <- assign_releves(
  x               = res,
  vegmatrix       = vm,
  mode            = "fuzzy",
  phi             = phi_cut,
  cover_transform = "sqrt",
  compare         = "max",
  min_cover       = 0,
  phi_mode        = "thresh",
  min_group_size  = 1
)
table(assign_fuzzy)
#> assign_fuzzy
#> c_1 c_3 
#>   2   2

# 9) Species for assigned (non-"not assigned") groups
labs_strict <- setdiff(names(table(assign_strict)), "not assigned")
if (length(labs_strict)) species_in_clusters(res, labs_strict)
#> $c_1
#> [1] "sp1" "sp5"
#> 
#> $c_3
#> [1] "sp2" "sp4" "sp6"

labs_fuzzy <- setdiff(names(table(assign_fuzzy)), "not assigned")
if (length(labs_fuzzy)) species_in_clusters(res_fuzzy, labs_fuzzy, min_phi = 0.15, top_k = 20)
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

------------------------------------------------------------------------

### (Optional) Attach assignments to header data frame

If you have a metadata table `hea` with a column `releve_number`, you
can join assignments back to it.  
(Example not evaluated in README builds.)

``` r
library(dplyr)
library(tibble)

assign_df <- tibble(
  plot_id         = names(assign_strict),
  cocktail_strict = unname(assign_strict),
  cocktail_fuzzy  = unname(assign_fuzzy)
)

hea2 <- hea %>%
  mutate(plot_id = as.character(releve_number)) %>%
  left_join(assign_df, by = "plot_id") %>%
  select(-plot_id)
```

------------------------------------------------------------------------

## Reference

See function help for details:

- `?cocktail_cluster` – deterministic Cocktail clustering  
- `?plot_cocktail` – draw dendrograms  
- `?clusters_at_cut` – find parent clusters at φ cuts  
- `?cocktail_fuzzy` – compute fuzzy species × node φ  
- `?assign_releves` – assign plots (strict/fuzzy)  
- `?species_in_clusters` – list diagnostic species for clusters

------------------------------------------------------------------------

© 2025 Denys Vynokurov & Helge Bruelheide. Licensed under MIT.
