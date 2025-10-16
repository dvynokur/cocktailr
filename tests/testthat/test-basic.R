test_that("cocktail_cluster works on a tiny matrix", {
  vm <- matrix(c(
    1,0,  # plot1: sp1 present
    0,1,  # plot2: sp2 present
    1,1   # plot3: both present
  ), nrow = 3, byrow = TRUE,
  dimnames = list(paste0("plot", 1:3), c("sp1", "sp2")))

  x <- cocktail_cluster(vm, progress = FALSE)

  expect_equal(ncol(x$Cluster.species), 2L)
  expect_equal(nrow(x$Cluster.species), 1L)   # n-1 merges

  labs <- clusters_at_cut(x, phi = 0.3)
  expect_true(is.character(labs))
})
