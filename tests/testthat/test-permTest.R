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
  mdl0 <- fitModel(
    spe = speSub,
    marks = "cellType",
    formula = as.formula(
        "Keratin_Tumour ~ 1"
    )
  )

  res <- permTest(mdl, mdl0, nsim = 10)

  expect_true(dplyr::between(res$pVal, 0, 1))
})
