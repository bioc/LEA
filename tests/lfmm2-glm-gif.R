# Standalone regression for column-wise GLM genomic control. Base R only;
# no optimized LFMM2 implementation or development test helpers are required.
library(LEA)

local({
    checks <- 0L
    check <- function(condition, message) {
        checks <<- checks + 1L
        if (!isTRUE(condition)) stop(message, call. = FALSE)
    }
    close <- function(actual, expected, message) {
        check(length(actual) == length(expected), paste(message, "length"))
        a <- as.numeric(actual)
        b <- as.numeric(expected)
        check(all(is.finite(a)) && all(is.finite(b)), paste(message, "finite"))
        check(all(abs(a - b) <= 1e-10 + 1e-7 * abs(b)), message)
    }

    # Match the established binary fixture, with p != d so the old per-locus
    # GIF vector cannot accidentally satisfy the covariate-count assertion.
    set.seed(7103L, kind = "Mersenne-Twister", normal.kind = "Inversion")
    n <- 80L
    p <- 40L
    K <- 2L
    X <- matrix(rnorm(n * 3L), n, 3L)
    A <- matrix(seq(0.15, 0.45, length.out = 3L * K), 3L, K)
    U <- X %*% A + matrix(rnorm(n * K), n, K)
    V <- sweep(matrix(rnorm(p * K), p, K), 2L,
               seq(1.5, 0.7, length.out = K), "*")
    B <- matrix(0, p, 3L)
    B[1:6, ] <- matrix(seq(-0.8, 0.9, length.out = 18L), 6L, 3L)
    probability <- plogis(0.25 * (X %*% t(B) + tcrossprod(U, V)))
    Y <- matrix(rbinom(n * p, 1L, probability), n, p)
    colnames(X) <- paste0("env", 1:3)

    for (d in c(3L, 1L)) {
        env <- X[, seq_len(d), drop = FALSE]
        fit <- lfmm2(Y, env, K)
        design <- data.frame(env, fit@U)
        for (family in list(binomial(), gaussian())) {
            label <- paste("d =", d, family$family)
            raw <- lfmm2.test(fit, Y, env, linear = FALSE,
                              family = family, genomic.control = FALSE)
            calibrated <- lfmm2.test(fit, Y, env, linear = FALSE,
                                     family = family, genomic.control = TRUE)

            # Derive statistics and uncalibrated inference independently using
            # ordinary per-locus GLMs, rather than another LEA calibration path.
            z <- probability <- matrix(NA_real_, p, d)
            for (j in seq_len(p)) {
                model <- glm(Y[, j] ~ ., data = design, family = family)
                coefficients <- summary(model)$coefficients
                z[j, ] <- coefficients[seq_len(d) + 1L, 3L]
                probability[j, ] <- coefficients[seq_len(d) + 1L, 4L]
            }
            close(raw$zscores, z, paste(label, "raw statistics"))
            close(raw$pvalues, probability, paste(label, "raw p-values"))
            check(identical(dim(raw$zscores), c(p, d)), paste(label, "raw shape"))
            check(is.null(raw$gif), paste(label, "uncalibrated GIF"))

            gif <- vapply(seq_len(d), function(k)
                median(z[, k]^2) / qchisq(0.5, df = 1L), numeric(1L))
            check(length(calibrated$gif) == d,
                  paste("GLM GIF count: expected", d, "got", length(calibrated$gif)))
            close(calibrated$gif, gif, paste(label, "column-wise GIF"))
            expected <- pchisq(sweep(z^2, 2L, gif, "/"), df = 1L,
                               lower.tail = FALSE)
            close(calibrated$pvalues, expected, paste(label, "calibrated p-values"))
            for (k in seq_len(d)) {
                close(calibrated$pvalues[, k],
                      pchisq(z[, k]^2 / gif[k], df = 1L, lower.tail = FALSE),
                      paste(label, "covariate", k))
            }
            check(identical(calibrated$zscores, raw$zscores),
                  paste(label, "statistics changed by calibration"))
            check(identical(dim(calibrated$pvalues), c(p, d)),
                  paste(label, "calibrated p-value shape"))
            check(is.null(calibrated$fscores) && is.null(calibrated$adj.r.squared),
                  paste(label, "GLM result fields"))
            check(identical(names(calibrated), names(raw)),
                  paste(label, "result names"))
            if (d == 1L) {
                # Historical single-covariate behavior uses one pooled median;
                # it must remain unchanged rather than use a new shape rule.
                historical <- median(raw$zscores^2) / qchisq(0.5, df = 1L)
                check(identical(calibrated$gif, historical),
                      paste(label, "single-covariate GIF changed"))
                check(identical(calibrated$pvalues,
                                pchisq(raw$zscores^2 / historical, df = 1L,
                                       lower.tail = FALSE)),
                      paste(label, "single-covariate p-values changed"))
            }
        }
    }
    cat("LFMM2 GLM GIF regression:", checks,
        "checks passed (d=3/d=1; binomial/Gaussian).\n")
})
