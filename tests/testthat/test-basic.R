test_that("cocktail runs on dune", {
  data(dune, package = "vegan")
  x <- cocktail_cluster(dune)
  expect_equal(ncol(x$Cluster.species), ncol(dune))
  expect_equal(nrow(x$Cluster.species), ncol(dune) - 1)
  cg <- cut_groups(x, phi = 0.3)
  expect_true(all(c("species","group") %in% names(cg)))
})
