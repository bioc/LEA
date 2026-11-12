test_that("linear t statistics and p-values agree with per-locus lm", {
    for (dat in lfmm_fixtures[c("simple", "multi")]) {
        fit <- lfmm2(dat$Y, dat$X, dat$K)
        ans <- lfmm2.test(fit, dat$Y, dat$X, genomic.control = FALSE)
        ref <- reference_association(dat$Y, dat$X, fit@U)
        expect_numeric_close(ans$zscores, ref$stat)
        expect_p_close(ans$pvalues, ref$p)
        expect_result_close(ans, historical_test(fit, dat$Y, dat$X,
                                                 genomic.control = FALSE))
        expect_identical(names(ans), c("pvalues", "zscores", "fscores",
                                       "adj.r.squared", "gif"))
        expect_null(ans$fscores)
        expect_null(ans$adj.r.squared)
        expect_null(ans$gif)
    }
})

test_that("linear historical dimensions and response names are retained", {
    for (dat in lfmm_fixtures[c("simple", "multi")]) {
        colnames(dat$Y) <- paste0("locus", seq_len(ncol(dat$Y)))
        colnames(dat$X) <- paste0("env", seq_len(ncol(dat$X)))
        fit <- lfmm2(dat$Y, dat$X, dat$K)
        raw <- lfmm2.test(fit, dat$Y, dat$X, genomic.control = FALSE)
        cal <- lfmm2.test(fit, dat$Y, dat$X)
        p <- ncol(dat$Y)
        d <- ncol(dat$X)
        if (d == 1L) {
            expect_null(dim(raw$pvalues))
            expect_identical(dim(raw$zscores), c(p, 1L))
            expect_identical(dim(cal$pvalues), c(p, 1L))
            expect_identical(names(raw$pvalues), paste("Response", colnames(dat$Y)))
        } else {
            expect_identical(dim(raw$pvalues), c(d, p))
            expect_identical(dim(raw$zscores), c(d, p))
            expect_identical(rownames(raw$zscores), colnames(dat$X))
            expect_identical(colnames(raw$zscores), paste("Response", colnames(dat$Y)))
            expect_identical(names(cal$gif), colnames(dat$X))
        }
        expect_result_close(cal, historical_test(fit, dat$Y, dat$X))
    }
})

test_that("linear tests support new loci without re-estimating factors", {
    dat <- lfmm_fixtures$multi
    fit <- lfmm2(dat$Y[, 1:40], dat$X, dat$K)
    Ynew <- dat$Y[, 41:55]
    ans <- lfmm2.test(fit, Ynew, dat$X, genomic.control = FALSE)
    ref <- reference_association(Ynew, dat$X, fit@U)
    expect_identical(dim(ans$zscores), c(3L, 15L))
    expect_numeric_close(ans$zscores, ref$stat)
    expect_p_close(ans$pvalues, ref$p)
    altered <- fit
    altered@B <- matrix(123)
    altered@V <- matrix(456)
    altered@K <- 99L
    altered@lambda <- 12
    expect_identical(lfmm2.test(altered, Ynew, dat$X, genomic.control = FALSE), ans)
})

test_that("association statistics respect paired row and locus permutations", {
    dat <- lfmm_fixtures$multi
    fit <- lfmm2(dat$Y, dat$X, dat$K)
    raw <- lfmm2.test(fit, dat$Y, dat$X, genomic.control = FALSE)
    rows <- rev(seq_len(nrow(dat$Y)))
    loci <- rev(seq_len(ncol(dat$Y)))
    perm <- fit
    perm@U <- fit@U[rows, , drop = FALSE]
    ans <- lfmm2.test(perm, dat$Y[rows, loci], dat$X[rows, ], genomic.control = FALSE)
    expect_numeric_close(ans$zscores, raw$zscores[, loci])
    expect_p_close(ans$pvalues, raw$pvalues[, loci])
})

test_that("p-value comparator checks extreme tails and exact underflow", {
    # A deterministic strong association actually exercises the logarithmic path.
    x <- seq(-1, 1, length.out = 60)
    U <- matrix(sin(seq_len(60)), ncol = 1)
    Y <- cbind(4 * x + 0.1 * cos(seq_len(60)),
               -3 * x + 0.15 * cos(2 * seq_len(60)))
    fit <- lfmm2(lfmm_fixtures$simple$Y, matrix(x), 1L)
    fit@U <- U
    ans <- lfmm2.test(fit, Y, matrix(x), genomic.control = FALSE)
    ref <- reference_association(Y, matrix(x), U)
    expect_true(all(ref$p < 1e-8 & ref$p > 0))
    expect_p_close(ans$pvalues, ref$p)
    expect_p_close(c(0, 1e-100), c(0, 1e-100))
})
