library("spatstat.explore")
spe <- SpatialDatasets::spe_Keren_2018()

test_that("summary.multipointRppm works", {
    speSub <- subset(spe, , imageID == "5")

    mdl <- fitModel(
        spe = speSub,
        model = "kppm",
        marks = "cellType",
        formula = as.formula("Keratin_Tumour ~ 1"),
        interaction = "Hardcore",
        threshold = 10
    )
    expect_equal(is(mdl), "multipointRppm")
    ret1 <- summary(mdl)
    ret2 <- spatstat.model::summary.kppm(mdl)
    expect_false(identical(ret1, ret2))
    #test written by claude Opus 5.5
    expect_true(all.equal(ret1$coefs.SE.CI$S.E.,
          unname(sqrt(diag(vcov(as.ppm(mdl) |> (\(x) mdl)(), fast = TRUE))))))
})

test_that("summary.multipointRppm doesn't change anything if not kppm", {
    speSub <- subset(spe, , imageID == "5")

    mdl <- fitModel(
        spe = speSub,
        model = "ppm",
        marks = "cellType",
        formula = as.formula("Keratin_Tumour ~ 1"),
        interaction = "Hardcore",
        threshold = 10
    )
    expect_equal(is(mdl), "multipointRppm")
    ret1 <- summary(mdl)
    ret2 <- spatstat.model::summary.ppm(mdl)
    expect_true(identical(ret1, ret2))
})

test_that("summary.multipointRppm doesn't change anything if not kppm", {
    speSub <- subset(spe, , imageID == "5")

    mdl <- fitModel(
        spe = speSub,
        model = "kppm",
        marks = "cellType",
        formula = as.formula("Keratin_Tumour ~ 1"),
        interaction = "Hardcore",
        threshold = 10
    )
    expect_equal(is(mdl), "multipointRppm")
    ret1 <- summary(mdl, fast = FALSE)
    ret2 <- spatstat.model::summary.kppm(mdl)
    expect_true(identical(ret1, ret2))
})