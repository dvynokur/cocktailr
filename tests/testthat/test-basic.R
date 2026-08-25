test_that("cocktailr core workflow works on a tiny matrix", {
  vm <- matrix(c(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    1, 1, 0,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  ), nrow = 7, byrow = TRUE,
  dimnames = list(paste0("plot", 1:7), c("sp1", "sp2", "sp3")))

  x <- cocktail_cluster(vm, progress = FALSE)

  expect_equal(ncol(x$Cluster.species), 3L)
  expect_equal(nrow(x$Cluster.species), 2L)   # n - 1 merges

  labs <- clusters_at_cut(x, phi_cut = 0.3)
  expect_true(is.character(labs))

  sel <- select_clusters(
    x,
    min_phi = -1,
    min_k = 1,
    min_n = 1,
    min_score = -Inf,
    score_method = "h",
    return = "table"
  )

  expect_s3_class(sel, "data.frame")
  expect_true(all(c("cluster", "id", "h", "k", "m", "n", "score") %in% names(sel)))

  d_cont <- cluster_dist(
    x,
    method = "containment",
    return = "dist"
  )

  d_phi <- cluster_dist(
    x,
    method = "phi",
    return = "dist"
  )

  expect_s3_class(d_cont, "dist")
  expect_s3_class(d_phi, "dist")

  tab <- cluster_dist(
    x,
    method = "containment",
    return = "table"
  )

  expect_s3_class(tab, "data.frame")
  expect_true(all(c(
    "cluster_1", "cluster_2",
    "id_1", "id_2",
    "n_1", "n_2",
    "shared_n",
    "containment_1_in_2",
    "containment_2_in_1",
    "containment_max",
    "containment_min",
    "jaccard",
    "phi",
    "method",
    "distance"
  ) %in% names(tab)))

  expect_true(all(tab$method == "containment"))
  expect_true(all(tab$distance >= 0))
})
