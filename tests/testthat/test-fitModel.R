library("spatstat.explore")
spe <- SpatialDatasets::spe_Keren_2018()

test_that("fitModel works with treshold", {
  speSub <- subset(spe, , imageID == "5")

  mdl <- fitModel(spe = speSub,
                  marks = "cellType",
                  formula = as.formula("Keratin_Tumour ~ 1"),
                  threshold = 10)
  expect_equal(is(mdl), "ppm")
})

test_that("fitModel fails with incorrect formula", {
  speSub <- subset(spe, , imageID == "5")

  expect_error(fitModel(spe = speSub,
                  marks = "cellType",
                  formula = as.formula("Keratin_Tumour"),
                  threshold = 10))
})

test_that("fitModel works with mixed formulas of changed and unchanged functions", {
  speSub <- subset(spe, , imageID == "5")

  mdl <- fitModel(spe = speSub,
                  marks = "cellType",
                  formula = as.formula("Keratin_Tumour ~ s(x) + log(density.ppp(Endothelial))"),
                  threshold = 10)
  expect_equal(is(mdl), "ppm")
})