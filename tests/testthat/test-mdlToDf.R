spe <- SpatialDatasets::spe_Keren_2018()

test_that("mdlToDf works", {
    mdlLs <- fitModelAcrossImages(
        spe = spe,
        imageId = "imageID",
        imageLs = list("1", "2"),
        marks = "cellType",
        sharedModel = FALSE,
        formula = as.formula("Keratin_Tumour ~
          spatstat.geom::distfun(CD8_T_cell)"),
        threshold = 10
    )
    imageCovariates <- c("imageID", "tumour_type")
    mdlDf <- mdlToDf(
        mdlLs = mdlLs,
        imageCovariates = imageCovariates
    )
    expect_true(nrow(mdlDf) == length(unique(mdlDf$imageID)) *
        (length(imageCovariates) + 1))
})
