# Compare values AND ordered conditions, including cases deliberately retained
# on the historical summary path. The frozen oracle remains unchanged.
capture_lfmm_conditions <- function(expr) {
    warnings <- character()
    error <- NULL
    value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
    list(value = value, warnings = warnings, error = error)
}

expect_lfmm_compatible <- function(fit, Y, X, full = FALSE) {
    for (control in c(FALSE, TRUE)) {
        actual <- capture_lfmm_conditions(lfmm2.test(fit, Y, X, full = full,
                                                    genomic.control = control))
        old <- capture_lfmm_conditions(historical_test(fit, Y, X, full = full,
                                                        genomic.control = control))
        expect_identical(actual$error, old$error)
        expect_identical(actual$warnings, old$warnings)
        if (is.null(old$error)) expect_result_close(actual$value, old$value)
    }
}

optimization_edge_cases <- function() {
    dat <- lfmm_fixtures$multi
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    Y <- dat$Y[, 1:6]
    named <- Y
    colnames(named) <- c("", "duplicate", "duplicate", "locus 4", "x", "y")
    alias_U <- fit
    alias_U@U[, 2] <- alias_U@U[, 1]
    missing_U <- fit
    missing_U@U[1, 1] <- NA
    small <- fit
    small@U <- fit@U[1:6, , drop = FALSE]
    perfect <- cbind(dat$X, fit@U) %*% matrix(seq_len(30), ncol = 6)
    list(
        named = list(fit = fit, Y = named, X = dat$X),
        constants = list(fit = fit, Y = cbind(Y, zero = 0, constant = 7), X = dat$X),
        alias_X = list(fit = fit, Y = Y, X = cbind(dat$X, 2 * dat$X[, 1])),
        near_alias_X = list(fit = fit, Y = Y,
                            X = cbind(dat$X, dat$X[, 1] + 1e-10 * dat$X[, 2])),
        alias_U = list(fit = alias_U, Y = Y, X = dat$X),
        missing_U = list(fit = missing_U, Y = Y, X = dat$X),
        saturated = list(fit = small, Y = Y[1:6, ], X = dat$X[1:6, ]),
        perfect = list(fit = fit, Y = perfect, X = dat$X))
}

test_that("optimized linear inference preserves exceptional values and conditions", {
    for (case in optimization_edge_cases())
        expect_lfmm_compatible(case$fit, case$Y, case$X)
})

test_that("optimized linear inference preserves file-based NA omission", {
    dat <- lfmm_fixtures$simple
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    Y <- dat$Y[, 1:4]
    Y[1, 1] <- NA
    path <- tempfile(fileext = ".lfmm")
    on.exit(unlink(path))
    write.table(Y, path, row.names = FALSE, col.names = FALSE)
    expect_lfmm_compatible(fit, path, dat$X)
})

test_that("optimized full inference preserves exceptional values and conditions", {
    for (case in optimization_edge_cases())
        expect_lfmm_compatible(case$fit, case$Y, case$X, full = TRUE)
})

test_that("optimized full inference preserves file-based NA omission failures", {
    dat <- lfmm_fixtures$simple
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    Y <- dat$Y[, 1:4]
    Y[1, 1] <- NA
    path <- tempfile(fileext = ".lfmm")
    on.exit(unlink(path))
    write.table(Y, path, row.names = FALSE, col.names = FALSE)
    expect_lfmm_compatible(fit, path, dat$X, full = TRUE)
})

test_that("full nuisance QR reuse preserves residuals and historical model ranks", {
    dat <- lfmm_fixtures$multi
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    nuisance_Y <- lm(dat$Y ~ ., data = data.frame(fit@U))
    nuisance_X <- lm(dat$X ~ ., data = data.frame(fit@U))
    expect_numeric_close(qr.resid(nuisance_Y$qr, dat$X), nuisance_X$residuals)
    for (X in list(dat$X, matrix(1, nrow(dat$X), 1),
                   cbind(dat$X, 2 * dat$X[, 1])))
        expect_lfmm_compatible(fit, dat$Y[, 1:6], X, full = TRUE)
})
