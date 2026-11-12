test_that("linear genomic control uses global per-covariate medians", {
    for (dat in lfmm_fixtures[c("simple", "multi")]) {
        fit <- lfmm2(dat$Y, dat$X, dat$K)
        raw <- lfmm2.test(fit, dat$Y, dat$X, genomic.control = FALSE)
        ans <- lfmm2.test(fit, dat$Y, dat$X)
        gif <- if (ncol(dat$X) == 1L) median(raw$zscores^2) else
            apply(raw$zscores^2, 1, median)
        gif <- gif / qchisq(0.5, 1, lower.tail = FALSE)
        expect_length(ans$gif, ncol(dat$X))
        expect_numeric_close(ans$gif, gif)
        expect_p_close(ans$pvalues, pchisq(raw$zscores^2 / gif, 1, lower.tail = FALSE))
        expect_identical(ans$zscores, raw$zscores)
        explicit <- lfmm2.test(fit, dat$Y, dat$X, full = FALSE,
                              genomic.control = TRUE, linear = TRUE,
                              family = binomial(link = "logit"))
        expect_identical(ans, explicit)
    }
})

test_that("single-covariate GLM genomic control pools across loci", {
    dat <- lfmm_fixtures$binary
    X <- dat$X[, 1, drop = FALSE]
    fit <- lfmm2(dat$Y, X, dat$K)
    raw <- lfmm2.test(fit, dat$Y, X, linear = FALSE, genomic.control = FALSE)
    ans <- lfmm2.test(fit, dat$Y, X, linear = FALSE)
    gif <- median(raw$zscores^2) / qchisq(0.5, 1)
    expect_length(ans$gif, 1L)
    expect_numeric_close(ans$gif, gif)
    expect_p_close(ans$pvalues, pchisq(raw$zscores^2 / gif, 1, lower.tail = FALSE))
    expect_identical(ans, historical_test(fit, dat$Y, X, linear = FALSE))
    expect_identical(dim(ans$pvalues), c(ncol(dat$Y), 1L))
    expect_identical(dim(ans$zscores), c(ncol(dat$Y), 1L))
})

test_that("HISTORICAL: multivariate GLM calibration uses one GIF per locus", {
    dat <- lfmm_fixtures$binary
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    raw <- lfmm2.test(fit, dat$Y, dat$X, linear = FALSE, genomic.control = FALSE)
    # Preserve the defective calculation only in the frozen oracle. Production
    # calibration intentionally differs; raw statistics and inference do not.
    ans <- historical_test(fit, dat$Y, dat$X, linear = FALSE)
    historical_gif <- apply(raw$zscores^2, 1, median) / qchisq(0.5, 1)
    expect_length(ans$gif, ncol(dat$Y))
    expect_numeric_close(ans$gif, historical_gif, rtol = 1e-6, atol = 1e-7)
    expect_p_close(ans$pvalues,
                   pchisq(raw$zscores^2 / historical_gif, 1, lower.tail = FALSE))
    corrected <- lfmm2.test(fit, dat$Y, dat$X, linear = FALSE)
    expect_identical(corrected$zscores, ans$zscores)
    expect_gt(max(abs(corrected$pvalues - ans$pvalues)), 1e-6)
    expect_identical(raw, historical_test(fit, dat$Y, dat$X,
                                         linear = FALSE, genomic.control = FALSE))
})

test_that("MATHEMATICAL: intended GLM calibration reference pools by covariate", {
    dat <- lfmm_fixtures$binary
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    ref <- reference_association(dat$Y, dat$X, fit@U, binomial())
    z <- t(ref$stat) # loci by covariates
    gif <- apply(z^2, 2, median) / qchisq(0.5, 1)
    intended <- pchisq(sweep(z^2, 2, gif, "/"), 1, lower.tail = FALSE)
    expect_length(gif, ncol(dat$X))
    for (j in seq_len(ncol(dat$X))) {
        expect_numeric_close(gif[j], median(z[, j]^2) / qchisq(0.5, 1))
        expect_p_close(intended[, j], pchisq(z[, j]^2 / gif[j], 1, lower.tail = FALSE))
    }
})

test_that("MATHEMATICAL: multivariate GLM API returns one GIF per covariate", {
    dat <- lfmm_fixtures$binary
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    ref <- reference_association(dat$Y, dat$X, fit@U, binomial())
    z <- t(ref$stat)
    gif <- apply(z^2, 2, median) / qchisq(0.5, 1)
    ans <- lfmm2.test(fit, dat$Y, dat$X, linear = FALSE)
    expect_length(ans$gif, ncol(dat$X))
    expect_numeric_close(ans$gif, gif, rtol = 1e-6, atol = 1e-7)
    expect_p_close(ans$pvalues, pchisq(sweep(z^2, 2, gif, "/"), 1, lower.tail = FALSE))
    expect_identical(dim(ans$pvalues), c(ncol(dat$Y), ncol(dat$X)))
    expect_identical(dim(ans$zscores), c(ncol(dat$Y), ncol(dat$X)))
    raw <- lfmm2.test(fit, dat$Y, dat$X, linear = FALSE, genomic.control = FALSE)
    expect_identical(ans$zscores, raw$zscores)
    expect_null(ans$fscores)
    expect_null(ans$adj.r.squared)
    for (j in seq_len(ncol(dat$X))) {
        expect_numeric_close(ans$gif[j], median(z[, j]^2) / qchisq(0.5, 1))
        expect_p_close(ans$pvalues[, j], pchisq(z[, j]^2 / gif[j], 1,
                                               lower.tail = FALSE))
    }
})
