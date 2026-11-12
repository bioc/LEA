test_that("packaged vignette data preserve fitting and full association results", {
    env <- new.env(parent = emptyenv())
    data("offset_example", package = "LEA", envir = env)
    dat <- env$offset_example
    fit <- lfmm2(dat$geno, dat$env, 2L, effect.sizes = TRUE)
    old <- historical_lfmm2(dat$geno, dat$env, 2L, effect.sizes = TRUE)
    expect_numeric_close(tcrossprod(fit@U, fit@V), tcrossprod(old@U, old@V))
    expect_numeric_close(fit@B, old@B)
    expect_result_close(lfmm2.test(fit, dat$geno, dat$env, full = TRUE),
                         historical_test(fit, dat$geno, dat$env, full = TRUE))
})

test_that("effect sizes reproduce downstream genetic gap and offset outputs", {
    dat <- lfmm_fixtures$multi
    prediction <- sweep(dat$X, 2, c(0.2, -0.1, 0.3), "+")
    loci <- 1:20
    # Both public functions average B over K; compare against the frozen fit.
    B <- (historical_lfmm2(dat$Y, dat$X, 1L, effect.sizes = TRUE)@B +
          historical_lfmm2(dat$Y, dat$X, 2L, effect.sizes = TRUE)@B) / 2
    B <- B[loci, , drop = FALSE]
    difference <- (dat$X - prediction) %*% t(B)
    eig <- eigen(cov(B))
    for (fun in list(genetic.gap, genetic.offset)) {
        ans <- fun(dat$Y, dat$X, pred.env = prediction,
                   K = c(1L, 2L), candidate.loci = loci)
        expect_numeric_close(ans$offset, rowMeans(difference^2))
        expect_numeric_close(ans$distance, rowMeans(abs(difference)))
        expect_numeric_close(ans$eigenvalues, eig$values)
        # Eigenvector signs are arbitrary; reconstruct the covariance instead.
        expect_numeric_close(ans$vectors %*% diag(ans$eigenvalues) %*% t(ans$vectors),
                             cov(B))
    }
})
