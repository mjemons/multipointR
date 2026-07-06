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

test_that("speToPP returns a ppp object with a continuous mark", {
  spe <- SpatialDatasets::spe_Keren_2018()
  speSub <- subset(spe, , imageID == "6")
  colData(speSub)$Na <- assay(speSub)["Na",] %>% t %>% as.matrix %>% 
  data.frame()
  pp <- speToPPP(speSub, mark = "Na", continuous = TRUE)

  expect_true(is(pp, "ppp"))
})

test_that("speToPP returns a ppp object with a custom window", {
  spe <- SpatialDatasets::spe_Keren_2018()
  speSub <- subset(spe, , imageID == "6")

  pp <- speToPPP(speSub, mark = "cellType", window = as.owin(c(0,1,0,1)))

  expect_true(is(pp, "ppp"))
})
