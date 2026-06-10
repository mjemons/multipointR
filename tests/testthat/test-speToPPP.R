test_that("speToPP returns a ppp object", {
  spe <- SpatialDatasets::spe_Keren_2018()
  speSub <- subset(spe, , imageID == "6")
  pp <- speToPPP(spe, mark = "cellType")

  expect_true(is(pp, "ppp"))
})

test_that("speToPP returns correct number of points", {
  spe <- SpatialDatasets::spe_Keren_2018()
  speSub <- subset(spe, , imageID == "6")
  pp <- speToPPP(speSub, mark = "cellType")

  expect_true(spatstat.geom::npoints(pp) == nrow(colData(speSub)))
})
