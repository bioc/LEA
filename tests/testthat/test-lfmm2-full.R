test_that("HISTORICAL: full-test F, adjusted R-squared and p-values are preserved", {
    for (dat in lfmm_fixtures[c("simple", "multi")]) {
        fit <- lfmm2(dat$Y, dat$X, dat$K)
        ans <- lfmm2.test(fit, dat$Y, dat$X, full = TRUE, genomic.control = FALSE)
        ref <- reference_full_historical(dat$Y, dat$X, fit@U)
        expect_numeric_close(ans$fscores, ref$f)
        expect_numeric_close(ans$adj.r.squared, ref$adj)
        expect_p_close(ans$pvalues, ref$p)
        expect_null(ans$zscores)
        expect_null(ans$gif)
        expect_null(dim(ans$pvalues))
        expect_length(ans$fscores, ncol(dat$Y))
        expect_result_close(ans, historical_test(fit, dat$Y, dat$X,
                              full = TRUE, genomic.control = FALSE))
    }
})

test_that("HISTORICAL: full-test genomic control retains nominal degrees of freedom", {
    dat <- lfmm_fixtures$multi
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    raw <- lfmm2.test(fit, dat$Y, dat$X, full = TRUE, genomic.control = FALSE)
    ans <- lfmm2.test(fit, dat$Y, dat$X, full = TRUE)
    d <- ncol(dat$X)
    df2 <- nrow(dat$Y) - d - 1
    gif <- median(raw$fscores) / qf(0.5, d, df2, lower.tail = FALSE)
    expect_numeric_close(ans$gif, gif)
    expect_p_close(ans$pvalues, pf(raw$fscores / gif, d, df2, lower.tail = FALSE))
    expect_identical(ans$fscores, raw$fscores)
    expect_identical(ans$adj.r.squared, raw$adj.r.squared)
    expect_result_close(ans, historical_test(fit, dat$Y, dat$X, full = TRUE))
})

test_that("MATHEMATICAL: conventional partial F reference agrees with nested lm", {
    dat <- lfmm_fixtures$multi
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    # Include a rank-deficient environmental design to test effective ranks.
    for (X in list(dat$X, cbind(dat$X, 2 * dat$X[, 1]))) {
        ref <- reference_partial_f(dat$Y, X, fit@U)
        expect_identical(ref$df1, 3L)
        expect_identical(ref$df2, 74L)
        nuisance <- data.frame(fit@U)
        complete <- data.frame(fit@U, X)
        for (j in c(1L, 7L, 40L)) {
            y <- dat$Y[, j]
            nested <- anova(lm(y ~ ., data = nuisance), lm(y ~ ., data = complete))
            expect_numeric_close(ref$f[j], nested$F[2])
            expect_p_close(ref$p[j], nested$`Pr(>F)`[2])
        }
    }
})

test_that("HISTORICAL versus conventional full F differs by residual degrees of freedom", {
    dat <- lfmm_fixtures$multi
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    current <- lfmm2.test(fit, dat$Y, dat$X, full = TRUE, genomic.control = FALSE)
    partial <- reference_partial_f(dat$Y, dat$X, fit@U)
    legacy_df2 <- nrow(dat$Y) - ncol(dat$X) - 1
    expect_numeric_close(partial$f, current$fscores * partial$df2 / legacy_df2)
    expect_gt(max(abs(partial$p - current$pvalues)), 1e-6)
    # A possible rank-aware recalibration is computed, without any claim about
    # improved calibration or power; choosing it is a future scientific decision.
    gif <- median(partial$f) / qf(0.5, partial$df1, partial$df2)
    calibrated <- pf(partial$f / gif, partial$df1, partial$df2, lower.tail = FALSE)
    expect_length(calibrated, ncol(dat$Y))
    expect_true(all(is.finite(calibrated) & calibrated >= 0 & calibrated <= 1))
})
