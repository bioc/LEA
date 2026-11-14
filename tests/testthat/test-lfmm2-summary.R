test_that("Cb is an independent uncentered coefficient second moment in both modes", {
    dat <- lfmm_fixtures$multi
    for (lambda in c(1e-5, 0.1, 10)) {
        retained <- lfmm2(dat$Y, dat$X, dat$K, lambda, TRUE)
        compact <- lfmm2(dat$Y, dat$X, dat$K, lambda, FALSE)
        # Frozen pre-change oracle, independent of the summary and metric code.
        old <- historical_lfmm2(dat$Y, dat$X, dat$K, lambda, TRUE)
        expected <- crossprod(old@B) / ncol(dat$Y)
        expect_numeric_close(attr(retained, "gea.metric")$Cb, expected)
        expect_numeric_close(attr(compact, "gea.metric")$Cb, expected)
        expect_identical(compact@B, matrix(NA))
        expect_numeric_close(retained@B, old@B)
        expect_identical(lfmm2.test(retained, dat$Y, dat$X),
                         lfmm2.test(compact, dat$Y, dat$X))
    }
})

metric_object <- function(Cb) {
    object <- new("lfmm2Class", K = 1L, lambda = 0.1,
                  U = matrix(0, 4, 1), V = matrix(0, 3, 1), B = matrix(NA))
    attr(object, "gea.metric") <- list(version = 1L, Cb = Cb)
    object
}

test_that("analytical spectra reconstruct metrics and respect the 99 percent threshold", {
    rotation <- matrix(c(1, 1, -1, 1), 2) / sqrt(2)
    for (values in list(c(99, 1), c(2, 2), c(4, 0), c(0, 0), 3)) {
        P <- if (length(values) == 2L) rotation else matrix(1)
        Cb <- P %*% diag(values, nrow = length(values)) %*% t(P)
        s <- summary(metric_object(Cb))$gea.spectrum
        expect_true(s$available)
        expect_equal(s$eigenvalues, values, tolerance = 1e-12)
        # Reconstruction handles arbitrary eigenvector signs and degenerate bases.
        expect_equal(s$vectors %*% diag(s$eigenvalues, length(values)) %*%
                     t(s$vectors), Cb, ignore_attr = TRUE, tolerance = 1e-12)
        expect_equal(crossprod(s$vectors), diag(length(values)),
                     ignore_attr = TRUE, tolerance = 1e-12)
        if (sum(values) == 0) {
            expect_true(all(is.na(s$proportion)))
            expect_true(is.na(s$axes99))
        } else {
            expect_equal(s$proportion, values / sum(values))
            expect_equal(s$cumulative, cumsum(values / sum(values)))
            expect_identical(s$axes99, which(cumsum(values / sum(values)) >= .99)[1L])
        }
    }
    for (scale in c(1e-200, 1e200)) {
        s <- summary(metric_object(diag(c(3, 1)) * scale))$gea.spectrum
        expect_equal(s$proportion, c(.75, .25))
    }
})

test_that("historical coefficients use L normalization and never covariance", {
    for (B in list(matrix(rep(c(2, -3), each = 3), 3, 2),
                   matrix(c(2, 3, 4, 1, 2, 3), 3, 2), matrix(7, 1, 1))) {
        object <- metric_object(matrix(0))
        attr(object, "gea.metric") <- NULL
        object@B <- B; object@V <- matrix(0, nrow(B), 1)
        s <- summary(object)$gea.spectrum
        expected <- Reduce(`+`, lapply(seq_len(nrow(B)), function(i)
            outer(B[i, ], B[i, ]))) / nrow(B)
        expect_equal(s$vectors %*% diag(s$eigenvalues, ncol(B)) %*% t(s$vectors),
                     expected, ignore_attr = TRUE, tolerance = 1e-12)
        if (nrow(B) > 1L) expect_false(isTRUE(all.equal(expected, cov(B))))
    }
    for (file in c("effect.rds", "no_effect.rds", "single.rds")) {
        old <- readRDS(test_path("..", "lfmm2-s4", file))
        expect_identical(summary(old)$gea.spectrum$available,
                         all(is.finite(old@B)))
        expect_output(show(old), "provenance not recorded")
    }
})

test_that("metadata captures original and transformed inputs with exact predictor ordering", {
    dat <- lfmm_fixtures$multi
    geno <- dat$Y
    climate <- data.frame(a = dat$X[, 1] + 5, b = dat$X[, 2], unused = 1)
    fit <- lfmm2(geno, climate, 2L, formula = ~ b + log(a))
    info <- attr(fit, "fit.info")
    expect_identical(info$genotype$label, "geno")
    expect_identical(info$environment$label, "climate")
    expect_identical(info$environment$dimensions, c(80L, 3L))
    expect_identical(info$fitted.dimensions, c(80L, 2L))
    expect_identical(info$predictors, c("b", "log(a)"))
    expect_equal(info$predictor.sd, c(b = sd(climate$b), "log(a)" = sd(log(climate$a))))
    expect_identical(rownames(attr(fit, "gea.metric")$Cb), info$predictors)
    for (env in list(as.matrix(climate), climate)) {
        object <- lfmm2(geno, env, 2L)
        expect_identical(attr(object, "fit.info")$environment$mode,
                         if (is.data.frame(env)) "data.frame" else "matrix")
    }
    gfile <- tempfile(fileext = ".lfmm"); efile <- tempfile(fileext = ".env")
    on.exit(unlink(c(gfile, efile)), add = TRUE)
    write.table(geno, gfile, row.names = FALSE, col.names = FALSE)
    write.table(as.matrix(climate), efile, row.names = FALSE, col.names = FALSE)
    object <- lfmm2(gfile, efile, 2L)
    expect_identical(attr(object, "fit.info")$genotype$label, gfile)
    expect_identical(attr(object, "fit.info")$environment$mode, "file")
    expect_identical(attr(object, "fit.info")$environment$dimensions, dim(climate))
    # Capture the expression without evaluating it twice.
    count <- 0L
    object <- lfmm2({count <- count + 1L; geno}, climate, 2L)
    expect_identical(count, 1L)
})

test_that("malformed optional attributes fail informatively without breaking display", {
    for (metric in list(list(version = 2L), list(version = 1L, Cb = NA),
        list(version = 1L, Cb = matrix(Inf)),
        list(version = 1L, Cb = diag(c(1, -1))),
        list(version = 1L, Cb = matrix(c(1, 0, 1, 1), 2)))) {
        object <- metric_object(matrix(0))
        attr(object, "gea.metric") <- metric
        attr(object, "fit.info") <- list(version = 1L, genotype = "bad")
        expect_false(summary(object)$gea.spectrum$available)
        expect_output(show(object), "LFMM2 fit")
    }
    object <- metric_object(diag(c(1, -1e-16)))
    expect_true(summary(object)$gea.spectrum$available)
    expect_equal(summary(object)$gea.spectrum$eigenvalues, c(1, 0))
})

test_that("display dispatch is concise and RDS round trips preserve summaries", {
    dat <- lfmm_fixtures$simple
    for (retain in c(FALSE, TRUE)) {
        object <- lfmm2(dat$Y, dat$X, 1L, effect.sizes = retain)
        expect_identical(capture.output(show(object)), capture.output(print(object)))
        expect_identical(withVisible(print(object))$visible, FALSE)
        s <- summary(object)
        expect_s3_class(s, "summary.lfmm2")
        expect_output(print(s), "sd = 1", fixed = TRUE)
        expect_true(s$gea.spectrum$available)
        path <- tempfile(); on.exit(unlink(path), add = TRUE)
        saveRDS(object, path)
        expect_identical(summary(readRDS(path)), s)
        expect_identical(slotNames(object), c("K", "lambda", "B", "U", "V"))
    }
    single <- lfmm2(dat$Y[, 1, drop = FALSE], dat$X, 1L)
    expect_identical(dim(attr(single, "gea.metric")$Cb), c(1L, 1L))
})

test_that("blocked accumulation agrees independently across locus boundaries", {
    dat <- make_lfmm_fixture(20L, 100003L, 2L, 1L, 1921L)
    object <- lfmm2(dat$Y, dat$X, 1L, 0.1, FALSE)
    old <- historical_lfmm2(dat$Y, dat$X, 1L, 0.1, TRUE)
    expect_numeric_close(attr(object, "gea.metric")$Cb,
                         crossprod(old@B) / ncol(dat$Y))
    expect_identical(object@B, matrix(NA))
})
