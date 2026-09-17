spe <- SpatialDatasets::spe_Keren_2018()
speSub <- subset(spe, , imageID == "5")

test_that("permutation test produces valid p-values", {
  mdl <- fitModel(
    spe = speSub,
    marks = "cellType",
    formula = as.formula(
        "Keratin_Tumour ~ spatstat.geom::distfun(CD8_T_cell)"
    )
  )

  res <- permTest(mdl, "spatstat.geom::distfun(CD8_T_cell)", nsim = 10)

  expect_true(dplyr::between(res$pVal, 0, 1))
})
