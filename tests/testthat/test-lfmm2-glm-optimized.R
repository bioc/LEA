# The existing frozen helper retains the original per-locus glm/rbind path.
# Its raw GLM inference is unchanged by the earlier genomic-control correction.
glm_optimization_reference <- function(fit, Y, X, family, control) {
    result <- historical_test(fit, Y, X, linear = FALSE,
                              genomic.control = FALSE, family = family)
    if (control) {
        z <- result$zscores
        if (ncol(as.matrix(X)) == 1L) {
            result$gif <- median(z^2) / qchisq(0.5, 1, lower.tail = FALSE)
            result$pvalues <- pchisq(z^2 / result$gif, 1, lower.tail = FALSE)
        } else {
            result$gif <- apply(z^2, 2, median) / qchisq(0.5, 1, lower.tail = FALSE)
            result$pvalues <- pchisq(sweep(z^2, 2, result$gif, "/"), 1, lower.tail = FALSE)
        }
    }
    result
}

glm_optimization_capture <- function(expr) {
    warnings <- character()
    error <- NULL
    value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
    list(value = value, warnings = warnings, error = error)
}

expect_glm_optimization_compatible <- function(fit, Y, X, family = binomial()) {
    for (control in c(FALSE, TRUE)) {
        actual <- glm_optimization_capture(lfmm2.test(fit, Y, X, linear = FALSE,
                                      family = family, genomic.control = control))
        reference <- glm_optimization_capture(glm_optimization_reference(
                                      fit, Y, X, family, control))
        expect_identical(actual$error, reference$error)
        expect_identical(actual$warnings, reference$warnings)
        if (is.null(reference$error)) {
            expect_result_close(actual$value, reference$value, glm = TRUE)
            # Same glm.fit and summary.glm operations should also be bitwise
            # identical on the same platform, including all public attributes.
            expect_identical(actual$value, reference$value)
        }
    }
}

glm_optimization_families <- function() {
    list(binomial(), binomial("probit"), binomial("cauchit"), binomial("cloglog"),
         gaussian(), gaussian("log"), poisson(), quasipoisson(),
         quasibinomial(), Gamma("log"), Gamma("inverse"), inverse.gaussian("log"))
}

glm_optimization_response <- function(family, Y, binary) {
    if (family$family %in% c("binomial", "quasibinomial")) return(binary)
    # 9 is LEA's missing-genotype sentinel: exercise valid count responses,
    # rather than letting API comparisons stop at the common input validator.
    if (family$family %in% c("poisson", "quasipoisson")) return(round(abs(Y) * 2) %% 7)
    if (family$family %in% c("Gamma", "inverse.gaussian") || family$link != "identity")
        return(abs(Y) + 1)
    Y
}

test_that("shared GLM design preserves standard families and links exactly", {
    dat <- lfmm_fixtures$binary
    continuous <- lfmm_fixtures$multi$Y[, 1:4]
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    for (family in glm_optimization_families()) {
        Y <- glm_optimization_response(family, continuous, dat$Y[, 1:4])
        expect_true(LEA:::.lfmm2_glm_reusable_family(family))
        expect_glm_optimization_compatible(fit, Y, dat$X, family)
    }
})

test_that("direct glm.fit retains coefficients, SEs, weights and iteration controls", {
    skip_if_not("singular.ok" %in% names(formals(stats::glm.fit)),
                "This R version uses the original glm() fallback, not direct glm.fit().")
    dat <- lfmm_fixtures$binary
    continuous <- lfmm_fixtures$multi$Y[, 1:3]
    X <- dat$X
    rownames(X) <- paste0("individual", seq_len(nrow(X)))
    fit <- lfmm2(dat$Y, X, dat$K)
    for (family in glm_optimization_families()) {
        Y <- glm_optimization_response(family, continuous, dat$Y[, 1:3])
        design_data <- data.frame(X, fit@U)
        j <- 1L
        template <- glm(Y[, j] ~ ., data = design_data, family = family)
        design <- model.matrix(template)
        for (j in 2:3) {
            expected <- glm(Y[, j] ~ ., data = design_data, family = family)
            actual <- LEA:::.lfmm2_glm_refit(design, Y[, j], family, template)
            expect_identical(actual$coefficients, expected$coefficients)
            expect_identical(summary(actual)$coefficients, summary(expected)$coefficients)
            expect_identical(actual$rank, expected$rank)
            expect_identical(actual$df.residual, expected$df.residual)
            expect_identical(actual$iter, expected$iter)
            expect_identical(actual$converged, expected$converged)
            expect_identical(actual$boundary, expected$boundary)
            expect_identical(actual$prior.weights, expected$prior.weights)
            expect_identical(actual$weights, expected$weights)
            expect_identical(actual$y, expected$y)
            expect_identical(actual$offset, expected$offset)
            expect_identical(actual$control, expected$control)
            expect_identical(actual$model, expected$model)
        }
    }
})

test_that("GLM preallocation preserves d=1, single loci and result names", {
    dat <- lfmm_fixtures$binary
    colnames(dat$X) <- c("env one", "duplicate", "duplicate")
    colnames(dat$Y) <- paste0("locus", seq_len(ncol(dat$Y)))
    for (d in c(1L, 3L)) {
        X <- dat$X[, seq_len(d), drop = FALSE]
        fit <- lfmm2(dat$Y[, 1:20], X, dat$K)
        for (loci in list(21L, 21:25)) {
            expect_glm_optimization_compatible(fit, dat$Y[, loci, drop = FALSE], X)
        }
    }
})

test_that("aliased, constant and saturated GLM fits retain historical extraction", {
    dat <- lfmm_fixtures$binary
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    Y <- dat$Y[, 1:4]
    expect_glm_optimization_compatible(fit, Y, cbind(dat$X, 2 * dat$X[, 1]))
    expect_glm_optimization_compatible(fit, Y,
        cbind(dat$X, dat$X[, 1] + 1e-12 * dat$X[, 2]))
    aliased <- fit
    aliased@U[, 2] <- aliased@U[, 1]
    expect_glm_optimization_compatible(aliased, Y, dat$X)
    expect_glm_optimization_compatible(fit, cbind(Y, zero = 0, one = 1), dat$X)
    saturated <- fit
    saturated@U <- fit@U[1:5, , drop = FALSE]
    expect_glm_optimization_compatible(saturated, Y[1:5, ], dat$X[1:5, ])
    expect_glm_optimization_compatible(saturated, Y[1:5, ], dat$X[1:5, ], gaussian())
})

test_that("categorical contrasts and logical inputs keep formula semantics", {
    dat <- lfmm_fixtures$binary
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    X <- cbind(group = rep(c("north", "south", "west", "north"), 20),
               treatment = rep(c("high", "low"), 40))
    previous <- options(contrasts = c("contr.sum", "contr.poly"))
    on.exit(options(previous))
    expect_glm_optimization_compatible(fit, dat$Y[, 1:4], X)
    expect_glm_optimization_compatible(fit, dat$Y[, 1:4], dat$X[, 1, drop = FALSE] > 0)
})

test_that("GLM missing rows, nonfinite data and na.action behavior are preserved", {
    dat <- lfmm_fixtures$binary
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    missing <- fit
    missing@U[1, 1] <- NA
    expect_glm_optimization_compatible(missing, dat$Y[, 1:4], dat$X)
    nonfinite <- fit
    nonfinite@U[1, 1] <- Inf
    expect_glm_optimization_compatible(nonfinite, dat$Y[, 1:4], dat$X)
    path <- tempfile(fileext = ".lfmm")
    on.exit(unlink(path))
    Y <- dat$Y[, 1:4]
    Y[2:3, 2] <- NA
    write.table(Y, path, row.names = FALSE, col.names = FALSE)
    expect_glm_optimization_compatible(fit, path, dat$X)
    previous <- options(na.action = "na.exclude")
    on.exit(options(previous), add = TRUE)
    expect_glm_optimization_compatible(fit, path, dat$X)
    options(na.action = "na.fail")
    expect_glm_optimization_compatible(fit, path, dat$X)
    options(previous)
    Y[2:3, 2] <- Inf
    write.table(Y, path, row.names = FALSE, col.names = FALSE)
    expect_glm_optimization_compatible(fit, path, dat$X)
})

test_that("separated GLMs retain non-convergence and boundary diagnostics", {
    dat <- lfmm_fixtures$binary
    X <- matrix(seq(-4, 4, length.out = 80), ncol = 1)
    fit <- lfmm2(dat$Y, X, dat$K)
    Y <- cbind(as.integer(X[, 1] > 0), as.integer(X[, 1] > 0.5),
               as.integer(X[, 1] > -0.5))
    expect_glm_optimization_compatible(fit, Y, X)
    actual <- glm_optimization_capture(lfmm2.test(fit, Y, X, linear = FALSE,
                                                 genomic.control = FALSE))
    expect_gt(length(actual$warnings), 0L)
})

test_that("custom na.action keeps response-specific model frames and call counts", {
    dat <- lfmm_fixtures$multi
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    Y <- dat$Y[, 1:4]
    previous <- options(na.action = function(frame)
        frame[frame[[1L]] > 0, , drop = FALSE])
    on.exit(options(previous))
    expect_glm_optimization_compatible(fit, Y, dat$X, gaussian())
    counter <- new.env(parent = emptyenv())
    counter$n <- 0L
    options(na.action = function(frame) { counter$n <- counter$n + 1L; frame })
    actual <- lfmm2.test(fit, Y, dat$X, linear = FALSE, family = gaussian(),
                         genomic.control = FALSE)
    expect_identical(counter$n, ncol(Y))
    counter$n <- 0L
    expected <- historical_test(fit, Y, dat$X, linear = FALSE, family = gaussian(),
                                genomic.control = FALSE)
    expect_identical(counter$n, ncol(Y))
    expect_identical(actual, expected)
    for (action in list(stats::na.omit, stats::na.exclude, stats::na.fail,
                        stats::na.pass, NULL)) {
        options(na.action = action)
        expect_glm_optimization_compatible(fit, Y, dat$X, gaussian())
    }
})

test_that("fractional binomial responses retain initialization warnings", {
    dat <- lfmm_fixtures$binary
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    Y <- dat$Y[, 1:4] / 2
    expect_glm_optimization_compatible(fit, Y, dat$X)
    actual <- glm_optimization_capture(lfmm2.test(fit, Y, dat$X, linear = FALSE,
                                                 genomic.control = FALSE))
    expect_true(any(grepl("non-integer", actual$warnings)))
})

test_that("family factories, names and custom callbacks retain fallback behavior", {
    dat <- lfmm_fixtures$binary
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    Y <- dat$Y[, 1:4]
    for (family in c("binomial", "gaussian", "poisson"))
        expect_glm_optimization_compatible(fit, Y, dat$X, family)
    counter <- new.env(parent = emptyenv())
    counter$n <- 0L
    factory <- function() { counter$n <- counter$n + 1L; binomial() }
    expect_false(LEA:::.lfmm2_glm_reusable_family(factory))
    actual <- lfmm2.test(fit, Y, dat$X, linear = FALSE, family = factory,
                         genomic.control = FALSE)
    expect_identical(counter$n, ncol(Y))
    counter$n <- 0L
    expected <- historical_test(fit, Y, dat$X, linear = FALSE, family = factory,
                                genomic.control = FALSE)
    expect_identical(counter$n, ncol(Y))
    expect_identical(actual, expected)
    custom <- binomial()
    original <- custom$linkinv
    custom$linkinv <- function(eta) original(eta)
    expect_false(LEA:::.lfmm2_glm_reusable_family(custom))
    expect_glm_optimization_compatible(fit, Y, dat$X, custom)
    initializer <- binomial()
    initializer$initialize <- expression({
        n <- rep.int(1, nobs)
        mustart <- rep.int(0.25, nobs)
    })
    expect_false(LEA:::.lfmm2_glm_reusable_family(initializer))
    expect_glm_optimization_compatible(fit, Y, dat$X, initializer)
    fixed <- gaussian()
    fixed$dispersion <- 2
    expect_glm_optimization_compatible(fit, lfmm_fixtures$multi$Y[, 1:4], dat$X, fixed)
    other <- quasi(link = "log", variance = "mu^2")
    expect_false(LEA:::.lfmm2_glm_reusable_family(other))
    expect_glm_optimization_compatible(fit, abs(lfmm_fixtures$multi$Y[, 1:4]) + 1,
                                      dat$X, other)
    unusual_name <- binomial()
    unusual_name$family <- factor("binomial")
    expect_false(LEA:::.lfmm2_glm_reusable_family(unusual_name))
    expect_glm_optimization_compatible(fit, Y, dat$X, unusual_name)
})
