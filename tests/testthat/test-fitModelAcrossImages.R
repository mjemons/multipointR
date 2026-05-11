library("spatstat.explore")
spe <- SpatialDatasets::spe_Keren_2018()

test_that("fitModelAcrossImages works as single models", {
  mdlLs <- fitModelAcrossImages(spe = spe,
                imageId = "imageID",
                imageLs = list("1", "2"),
                sharedModel = FALSE,
                marks = "cellType",
                formula = as.formula("Keratin_Tumour ~ distfun(CD8_T_cell)"),
                threshold = 10)
  expect_equal(is(mdlLs)[[1]], "list")
  expect_true(length(mdlLs)>1)
})

test_that("fitModelAcrossImages works as single models with interactions", {
  mdlLs <- fitModelAcrossImages(spe = spe,
                imageId = "imageID",
                imageLs = list("1", "2"),
                sharedModel = FALSE,
                marks = "cellType",
                formula = as.formula("Keratin_Tumour ~ distfun(CD8_T_cell)"),
                threshold = 10)
  expect_true(!is.null(mdlLs[[1]]$interaction))
})

test_that("fitModelAcrossImages works as shared model with interactions", {
  mdl <- fitModelAcrossImages(spe = spe,
                imageId = "imageID",
                imageLs = list("1", "2"),
                sharedModel = TRUE,
                marks = "cellType",
                formula = as.formula("Keratin_Tumour ~ distfun(CD8_T_cell)"),
                threshold = 10)
  expect_true(!is.null(mdlLs$Inter$interaction))
})

test_that("fitModelAcrossImages works as shared model with random effects", {
  mdl <- fitModelAcrossImages(spe = spe,
                imageId = "imageID",
                imageLs = list("1", "2"),
                sharedModel = TRUE,
                marks = "cellType",
                formula = as.formula("Keratin_Tumour ~ distfun(CD8_T_cell) + (1|sample_id)"),
                threshold = 10)
  expect_true(!is.null(mdl))
})