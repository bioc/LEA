test_that("retained tied singular values are compared by reconstruction", {
    n <- 20L
    basis <- qr.Q(qr(cbind(1, outer(seq_len(n), 1:7, function(i, j) sin(i * j)))))
    X <- basis[, 8, drop = FALSE]
    # The leading tie is entirely retained; individual vectors may rotate.
    Y <- sweep(basis[, 2:6], 2, c(8, 8, 3, 2, 1), "*") %*%
        t(diag(8)[, 1:5])
    fit <- lfmm2(Y, X, 2L, lambda = 1, effect.sizes = TRUE)
    ref <- reference_fit(Y, X, 2L, 1)
    expect_numeric_close(tcrossprod(fit@U, fit@V), ref$W)
    expect_numeric_close(fit@B, ref$B)
    expect_numeric_close(crossprod(fit@V), diag(2))
})

test_that("zero-rank response has zero fitted reconstruction and effects", {
    X <- lfmm_fixtures$simple$X
    Y <- matrix(0, nrow(X), 4)
    fit <- lfmm2(Y, X, 1L, effect.sizes = TRUE)
    expect_numeric_close(tcrossprod(fit@U, fit@V), Y)
    expect_numeric_close(fit@B, matrix(0, 4, 1))
})

test_that("near-collinear covariates and extreme ridge values preserve historical fits", {
    dat <- lfmm_fixtures$simple
    X <- cbind(dat$X, dat$X + 1e-7 * sin(seq_len(nrow(dat$X))))
    # An independent inverse-square-root reference would itself be ill-conditioned
    # here. Use the fixed historical algorithm, without relaxing main tolerances.
    for (lambda in c(1e-8, 1e4)) {
        fit <- lfmm2(dat$Y, X, 1L, lambda, TRUE)
        old <- historical_lfmm2(dat$Y, X, 1L, lambda, TRUE)
        expect_numeric_close(tcrossprod(fit@U, fit@V), tcrossprod(old@U, old@V))
        expect_numeric_close(fit@B, old@B)
        expect_true(all(is.finite(fit@B)))
    }
})

test_that("HISTORICAL: fractional K is accepted and stored after truncation", {
    dat <- lfmm_fixtures$simple
    fit <- lfmm2(dat$Y, dat$X, 2.5)
    expect_identical(fit@K, 2L)
    old <- historical_lfmm2(dat$Y, dat$X, 2.5)
    expect_numeric_close(tcrossprod(fit@U, fit@V), tcrossprod(old@U, old@V))
    expect_error(lfmm2(dat$Y, dat$X, 1.5), "non-conformable")
})

test_that("invalid rank and nonfinite input failures are characterized", {
    dat <- lfmm_fixtures$simple
    for (K in list(0L, -1L, 100L, NA_integer_, integer())) {
        expect_error(lfmm2(dat$Y, dat$X, K))
        expect_error(historical_lfmm2(dat$Y, dat$X, K))
    }
    for (lambda in list(NA_real_, NaN, Inf, numeric(), c(1, 2))) {
        expect_error(lfmm2(dat$Y, dat$X, 1L, lambda))
        expect_error(historical_lfmm2(dat$Y, dat$X, 1L, lambda))
    }
    for (value in c(Inf, -Inf)) {
        Y <- dat$Y
        Y[1, 1] <- value
        expect_error(lfmm2(Y, dat$X, 1L), "infinite|missing")
        X <- dat$X
        X[1, 1] <- value
        expect_error(lfmm2(dat$Y, X, 1L), "infinite|missing")
    }
    expect_error(lfmm2(dat$Y[, FALSE, drop = FALSE], dat$X, 1L))
})
