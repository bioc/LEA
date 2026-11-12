test_that("binomial GLM statistics and tails match independent glm", {
    dat <- lfmm_fixtures$binary
    for (d in c(1L, 3L)) {
        X <- dat$X[, seq_len(d), drop = FALSE]
        fit <- lfmm2(dat$Y, X, dat$K)
        raw <- lfmm2.test(fit, dat$Y, X, linear = FALSE, genomic.control = FALSE)
        ref <- reference_association(dat$Y, X, fit@U, binomial())
        expect_numeric_close(raw$zscores, t(ref$stat), rtol = 1e-6, atol = 1e-7)
        expect_p_close(raw$pvalues, t(ref$p))
        expect_identical(dim(raw$zscores), c(ncol(dat$Y), d))
        expect_identical(dim(raw$pvalues), c(ncol(dat$Y), d))
        expect_null(raw$gif)
        expect_null(raw$fscores)
        expect_null(raw$adj.r.squared)
        expect_result_close(raw, historical_test(fit, dat$Y, X,
                            linear = FALSE, genomic.control = FALSE), glm = TRUE)
        expect_identical(raw, historical_test(fit, dat$Y, X,
                            linear = FALSE, genomic.control = FALSE))
    }
})

test_that("non-binomial families preserve family-specific inference", {
    dat <- lfmm_fixtures$simple
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    ans <- lfmm2.test(fit, dat$Y, dat$X, linear = FALSE,
                      genomic.control = FALSE, family = gaussian())
    ref <- reference_association(dat$Y, dat$X, fit@U, gaussian())
    expect_numeric_close(ans$zscores, t(ref$stat), rtol = 1e-6, atol = 1e-7)
    expect_p_close(ans$pvalues, t(ref$p))
    linear <- lfmm2.test(fit, dat$Y, dat$X, genomic.control = FALSE)
    expect_numeric_close(ans$zscores, linear$zscores, rtol = 1e-6, atol = 1e-7)
    expect_p_close(ans$pvalues, linear$pvalues)
})

test_that("GLM tests support single and previously unfitted loci", {
    dat <- lfmm_fixtures$binary
    fit <- lfmm2(dat$Y[, 1:20], dat$X, dat$K)
    for (loci in list(21L, 21:25)) {
        Y <- dat$Y[, loci, drop = FALSE]
        ans <- lfmm2.test(fit, Y, dat$X, linear = FALSE, genomic.control = FALSE)
        ref <- reference_association(Y, dat$X, fit@U, binomial())
        expect_identical(dim(ans$zscores), c(length(loci), 3L))
        expect_numeric_close(ans$zscores, t(ref$stat), rtol = 1e-6, atol = 1e-7)
        expect_p_close(ans$pvalues, t(ref$p))
    }
})
