spe <- SpatialDatasets::spe_Keren_2018()

test_that("fitModel works with treshold", {
  speSub <- subset(spe, , imageID == "15")

  mdl <- fitModel(spe = speSub,
                  marks = "cellType",
                  response = "Keratin_Tumour",
                  distanceTo = "CD8_T_cell",
                  threshold = 10)
})
