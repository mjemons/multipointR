spe <- SpatialDatasets::spe_Keren_2018()

test_that("fitModel works with treshold", {
    speSub <- subset(spe, , imageID == "5")

    mdl <- fitModel(
        spe = speSub,
        marks = "cellType",
        formula = as.formula("Keratin_Tumour ~ 1"),
        interaction = "Strauss",
        threshold = 10
    )
    expect_equal(is(mdl), "multipointRppm")

    p <- plot(mdl)
    is(p, "ggplot2::ggplot")
})
