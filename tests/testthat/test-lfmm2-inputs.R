test_that("matrix validation preserves explicit errors", {
    dat <- lfmm_fixtures$simple
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    for (sentinel in c(NA_real_, NaN, 9, -9)) {
        Y <- dat$Y
        Y[1, 1] <- sentinel
        expect_error(lfmm2(Y, dat$X, 1L), "missing values")
        expect_error(lfmm2.test(fit, Y, dat$X), "missing values")
    }
    expect_error(lfmm2(NULL, dat$X, 1L), "NULL value")
    expect_error(lfmm2(dat$Y, NULL, 1L), "NULL value")
    expect_error(lfmm2.test(fit, NULL, dat$X), "NULL value")
    expect_error(lfmm2.test(fit, dat$Y, NULL), "NULL value")
    X <- dat$X
    X[1, 1] <- NA
    expect_error(lfmm2(dat$Y, X, 1L), "environmental matrix contains NA")
    expect_error(lfmm2.test(fit, dat$Y, X), "environmental matrix contains NA")
    expect_error(lfmm2(dat$Y[-1, ], dat$X, 1L), "Number of rows")
    expect_error(lfmm2.test(fit, dat$Y[-1, ], dat$X), "Number of rows")
    expect_error(lfmm2(dat$Y, dat$X, c(1L, 2L)), "Multiple values")
    for (lambda in c(0, -1))
        expect_error(lfmm2(dat$Y, dat$X, 1L, lambda), "must be positive")
    wide <- matrix(seq_len(6), 2, 3)
    expect_error(lfmm2(dat$Y[1:2, ], wide, 1L), "more columns")
    expect_error(lfmm2.test(fit, dat$Y[1:2, ], wide), "more columns")
    expect_error(lfmm2.test(fit, dat$Y, dat$X, full = TRUE, linear = FALSE),
                 "only for linear models")
})

test_that("data frames and numeric environmental vectors retain coercion", {
    dat <- lfmm_fixtures$simple
    fit <- lfmm2(dat$Y, dat$X, 1L)
    other <- lfmm2(as.data.frame(dat$Y), as.vector(dat$X), 1L)
    expect_numeric_close(tcrossprod(other@U, other@V), tcrossprod(fit@U, fit@V))
    frame <- as.data.frame(dat$Y)
    ans <- lfmm2.test(fit, frame, as.vector(dat$X))
    raw <- lfmm2.test(fit, dat$Y, dat$X)
    expect_numeric_close(ans$zscores, raw$zscores)
    expect_p_close(ans$pvalues, raw$pvalues)
    expect_identical(rownames(ans$zscores), paste("Response", names(frame)))
    expect_result_close(ans, historical_test(fit, frame, as.vector(dat$X)))
})

test_that("file inputs agree with matrices and preserve extension and sentinel checks", {
    dat <- lfmm_fixtures$simple
    # Integer responses avoid decimal serialization obscuring import equivalence.
    Y <- lfmm_fixtures$binary$Y[1:60, 1:12]
    colnames(Y) <- paste0("V", seq_len(ncol(Y)))
    X <- matrix(seq_len(60) %% 7, ncol = 1)
    yp <- tempfile(fileext = ".lfmm")
    xp <- tempfile(fileext = ".env")
    on.exit(unlink(c(yp, xp)))
    write.table(Y, yp, row.names = FALSE, col.names = FALSE)
    write.table(X, xp, row.names = FALSE, col.names = FALSE)
    fit <- lfmm2(Y, X, 1L)
    filefit <- lfmm2(yp, xp, 1L)
    expect_numeric_close(tcrossprod(filefit@U, filefit@V), tcrossprod(fit@U, fit@V))
    expect_warning(ans <- lfmm2.test(filefit, yp, xp), "Reading large input files")
    expect_result_close(ans, lfmm2.test(fit, Y, X))
    expect_error(lfmm2(xp, xp, 1L), "extension")
    for (sentinel in c(9, -9)) {
        Y[1, 1] <- sentinel
        write.table(Y, yp, row.names = FALSE, col.names = FALSE)
        expect_error(lfmm2(yp, xp, 1L), "missing data")
        expect_warning(expect_error(lfmm2.test(fit, yp, xp), "missing data"),
                       "Reading large input files")
    }
})

test_that("HISTORICAL: single-locus linear and full testing fail during summary extraction", {
    dat <- lfmm_fixtures$simple
    fit <- lfmm2(dat$Y, dat$X, 1L)
    for (full in c(FALSE, TRUE))
        expect_error(lfmm2.test(fit, dat$Y[, 1, drop = FALSE], dat$X,
                               full = full, genomic.control = FALSE),
                     "atomic vectors")
})

test_that("HISTORICAL: a constant locus contaminates genomic-control medians", {
    dat <- lfmm_fixtures$simple
    fit <- lfmm2(dat$Y, dat$X, 1L)
    Y <- cbind(dat$Y[, 1:5], constant = 0)
    expect_warning(raw <- lfmm2.test(fit, Y, dat$X, genomic.control = FALSE),
                   "non-variable genomic sites")
    expect_true(all(is.finite(raw$pvalues[1:5])))
    expect_true(is.na(raw$zscores[6, 1]))
    expect_warning(cal <- lfmm2.test(fit, Y, dat$X), "non-variable genomic sites")
    expect_true(is.na(cal$gif))
    expect_true(all(is.na(cal$pvalues)))
    expect_warning(full <- lfmm2.test(fit, Y, dat$X, full = TRUE),
                   "non-variable genomic sites")
    expect_true(all(is.na(full$pvalues)))
})

test_that("HISTORICAL: aliased environmental coefficients can be confused with factors", {
    dat <- lfmm_fixtures$simple
    X <- cbind(dat$X, 2 * dat$X)
    fit <- lfmm2(dat$Y, X, 1L, lambda = 0.1)
    ans <- lfmm2.test(fit, dat$Y, X, genomic.control = FALSE)
    expect_result_close(ans, historical_test(fit, dat$Y, X, genomic.control = FALSE))
    # summary.lm omits the aliased X2, so positional extraction selects U1 as
    # the purported second environmental statistic instead of returning NA.
    m <- lm(dat$Y[, 1] ~ ., data = data.frame(X, fit@U))
    expect_true(is.na(coef(m)[3]))
    expect_numeric_close(ans$zscores[2, 1], tail(summary(m)$coefficients[, 3], 1))
})
