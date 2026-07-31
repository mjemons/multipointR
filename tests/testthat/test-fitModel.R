library("spatstat.explore")
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
    expect_true(verifyclass(mdl, "ppm"))
})

test_that("fitModel fails with incorrect formula", {
    speSub <- subset(spe, , imageID == "5")

    expect_error(fitModel(
        spe = speSub,
        marks = "cellType",
        formula = as.formula("Keratin_Tumour"),
        interaction = "StraussHard",
        threshold = 10
    ))
})

test_that("fitModel works with mixed formulas of changed
 and unchanged functions", {
    speSub <- subset(spe, , imageID == "5")

    mdl <- fitModel(
        spe = speSub,
        marks = "cellType",
        formula = 
        as.formula("Keratin_Tumour ~ sqrt(x) + log(density.ppp(Endothelial))"),
        interaction = "Fiksel",
        cellspacing = 1,
        threshold = 10
    )
    expect_equal(is(mdl), "multipointRppm")
    expect_true(verifyclass(mdl, "ppm"))
})

test_that("fitModel returns NULL as model if the covariate is absent", {
    speSub <- subset(spe, , imageID == "15")

    mdl <- fitModel(
        spe = speSub,
        marks = "cellType",
        formula = 
        as.formula("Keratin_Tumour ~ sqrt(x) + log(density.ppp(CD8_T_cell))"),
        threshold = 10
    )
    expect_true(is.null(mdl))
})

test_that("fitModel works with interactions", {
    speSub <- subset(spe, , imageID == "15")

    mdl <- fitModel(
        spe = speSub,
        marks = "cellType",
        formula = 
        as.formula("Keratin_Tumour ~ sqrt(x) + log(density.ppp(Endothelial))"),
        interaction = "Hardcore",
        threshold = 10
    )
    expect_true(!is.null(mdl$interaction))
})

test_that("fitModel works with offset", {
    speSub <- subset(spe, , imageID == "15")

    mdl <- fitModel(
        spe = speSub,
        marks = "cellType",
        formula =
        as.formula("Keratin_Tumour ~ offset(log(lambda)) + distfun(Endothelial)"),
        interaction = NULL
    )
    expect_true(!is.null(mdl$trend))
})

test_that("fitModel fails if interaction is not implemented", {
    speSub <- subset(spe, , imageID == "15")

    expect_error(mdl <- fitModel(
        spe = speSub,
        marks = "cellType",
        formula = 
        as.formula("Keratin_Tumour ~ offset(log(lambda)) + distfun(Endothelial)"),
        interaction = "Dirichlet"
    ))
})

test_that("fitModel works with interaction", {
    speSub <- subset(spe, , imageID == "15")

    mdl <- fitModel(
        spe = speSub,
        marks = "cellType",
        formula = as.formula("Keratin_Tumour ~ sqrt(x):distfun(Endothelial)"),
        interaction = "StraussHard"
    )
    expect_equal(is(mdl), "multipointRppm")
    expect_true(verifyclass(mdl, "ppm"))
})

test_that("fitModel works with sf polygon from `sosta`", {
    speSub <- subset(spe, , imageID == "15")

    sfeSub <- SpatialFeatureExperiment::toSpatialFeatureExperiment(speSub)
    (struct <- sosta::reconstructShapeDensityImage(
        sfeSub,
        marks = "cellType",
        markSelect = c("Keratin_Tumour")
    ))

    SpatialFeatureExperiment::annotGeometry(sfeSub, "tumour_mask") <- struct

    mdl <- fitModel(
        spe = sfeSub,
        marks = "cellType",
        formula = as.formula("Endothelial ~ tumour_mask")
    )
    expect_equal(is(mdl), "multipointRppm")
    expect_true(verifyclass(mdl, "ppm"))
})

test_that("fitModel works with enet lasso reg", {
    speSub <- subset(spe, , imageID == "15")

    mdl <- fitModel(
        spe = speSub,
        marks = "cellType",
        formula = as.formula("Keratin_Tumour ~ log(lambda) + splines::bs(x)"), 
        improve.type = "enet",
        improve.args = list(alpha = 1),
        interaction = "Hardcore"
    )
    expect_equal(is(mdl), "multipointRppm")
    expect_true(verifyclass(mdl, "ppm"))
})

test_that("fitModel works with enet lasso reg and relaxed fitting", {
    speSub <- subset(spe, , imageID == "15")

    mdl <- fitModel(
        spe = speSub,
        marks = "cellType",
        formula = as.formula("Keratin_Tumour ~ log(lambda) + splines::bs(x)"), 
        improve.type = "enet",
        improve.args = list(alpha = 1),
        interaction = "Hardcore"
    )

    mdlRelax <- fitModel(
        spe = speSub,
        marks = "cellType",
        formula = as.formula("Keratin_Tumour ~ log(lambda) + splines::bs(x)"), 
        improve.type = "enet",
        improve.args = list(alpha = 1),
        interaction = "Hardcore",
        relaxed = TRUE
    )
    expect_equal(is(mdlRelax), "multipointRppm")
    expect_true(verifyclass(mdlRelax, "ppm"))
    expect_true(length(coef(mdlRelax)) < length(coef(mdl)))
})
