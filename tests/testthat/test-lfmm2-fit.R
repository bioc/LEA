test_that("fitting preserves S4 slots, defaults and dimensions", {
    for (dat in lfmm_fixtures[c("simple", "multi")]) {
        fit <- lfmm2(dat$Y, dat$X, dat$K)
        expect_s4_class(fit, "lfmm2Class")
        expect_identical(slotNames(fit), c("K", "lambda", "B", "U", "V"))
        expect_identical(fit@K, as.integer(dat$K))
        expect_identical(fit@lambda, 1e-5)
        expect_identical(fit@B, matrix(NA))
        expect_identical(dim(fit@U), c(nrow(dat$Y), dat$K))
        expect_identical(dim(fit@V), c(ncol(dat$Y), dat$K))
        expect_true(validObject(fit))
        explicit <- lfmm2(dat$Y, dat$X, dat$K, lambda = 1e-5,
                          effect.sizes = FALSE)
        expect_equal(fit, explicit, tolerance = 0)
    }
})

test_that("historical reconstruction and B survive K and ridge choices", {
    for (dat in lfmm_fixtures[c("simple", "multi")]) {
        for (K in c(1L, 2L)) for (lambda in c(1e-5, 0.1, 10)) {
            fit <- lfmm2(dat$Y, dat$X, K, lambda, effect.sizes = TRUE)
            old <- historical_lfmm2(dat$Y, dat$X, K, lambda, TRUE)
            expect_identical(dim(fit@B), c(ncol(dat$Y), ncol(dat$X)))
            expect_identical(fit@lambda, lambda)
            expect_numeric_close(tcrossprod(fit@U, fit@V), tcrossprod(old@U, old@V))
            expect_numeric_close(fit@B, old@B)
            expect_numeric_close(crossprod(fit@V), diag(K))
            expect_numeric_close(tcrossprod(fit@U), tcrossprod(old@U))
            no_effect <- lfmm2(dat$Y, dat$X, K, lambda, FALSE)
            expect_identical(no_effect@B, matrix(NA))
            expect_numeric_close(tcrossprod(no_effect@U, no_effect@V),
                                 tcrossprod(fit@U, fit@V))
        }
    }
})

test_that("symmetric ridge reference independently reproduces W and B", {
    for (dat in lfmm_fixtures[c("simple", "multi")]) {
        for (lambda in c(0.1, 10)) {
            fit <- lfmm2(dat$Y, dat$X, dat$K, lambda, TRUE)
            ref <- reference_fit(dat$Y, dat$X, dat$K, lambda)
            expect_numeric_close(tcrossprod(fit@U, fit@V), ref$W)
            expect_numeric_close(fit@B, ref$B)
        }
    }
})

test_that("centering and row, locus and covariate permutations preserve fits", {
    dat <- lfmm_fixtures$multi
    fit <- lfmm2(dat$Y, dat$X, dat$K, lambda = 0.1, effect.sizes = TRUE)
    W <- tcrossprod(fit@U, fit@V)
    shifted <- lfmm2(sweep(dat$Y, 2, seq_len(ncol(dat$Y)) / 10, "+"),
                     sweep(dat$X, 2, c(2, -3, 4), "+"), dat$K, 0.1, TRUE)
    expect_numeric_close(tcrossprod(shifted@U, shifted@V), W)
    expect_numeric_close(shifted@B, fit@B)
    rows <- rev(seq_len(nrow(dat$Y)))
    loci <- rev(seq_len(ncol(dat$Y)))
    covariates <- c(3L, 1L, 2L)
    perm <- lfmm2(dat$Y[rows, loci], dat$X[rows, covariates], dat$K, 0.1, TRUE)
    expect_numeric_close(tcrossprod(perm@U, perm@V), W[rows, loci])
    expect_numeric_close(perm@B, fit@B[loci, covariates])
    # An orthogonal environmental rotation preserves the isotropic ridge fit.
    rotation <- qr.Q(qr(matrix(c(1, 2, 0, -2, 1, 1, 1, 0, 2), 3)))
    rotated <- lfmm2(dat$Y, dat$X %*% rotation, dat$K, 0.1, TRUE)
    expect_numeric_close(tcrossprod(rotated@U, rotated@V), W)
    expect_numeric_close(rotated@B, fit@B %*% rotation)
})

test_that("single-locus, constant and collinear fitting are characterized", {
    dat <- lfmm_fixtures$simple
    single <- lfmm2(dat$Y[, 1, drop = FALSE], dat$X, 1L, 0.1, TRUE)
    expect_identical(dim(single@V), c(1L, 1L))
    expect_identical(dim(single@B), c(1L, 1L))
    expect_numeric_close(tcrossprod(single@U, single@V),
                         scale(dat$Y[, 1, drop = FALSE], scale = FALSE))
    Y <- cbind(dat$Y, constant = 0)
    X <- cbind(dat$X, 2 * dat$X)
    fit <- lfmm2(Y, X, 1L, 0.1, TRUE)
    ref <- reference_fit(Y, X, 1L, 0.1)
    expect_numeric_close(tcrossprod(fit@U, fit@V), ref$W)
    expect_numeric_close(fit@B, ref$B)
    expect_numeric_close(fit@B[nrow(fit@B), ], c(0, 0))
})

test_that("fixtures are reproducible and numerical calls do not mutate inputs or RNG", {
    dat <- lfmm_fixtures$simple
    # Use a serialized copy so a future native in-place mutation cannot alter
    # both sides of the assertion through shared R storage.
    original <- unserialize(serialize(dat, NULL))
    set.seed(42)
    before <- .Random.seed
    recreated <- make_lfmm_fixture(60L, 40L, 1L, 1L, 7101L)
    expect_identical(.Random.seed, before)
    expect_identical(recreated, dat)
    fit <- lfmm2(dat$Y, dat$X, dat$K, effect.sizes = TRUE)
    lfmm2.test(fit, dat$Y, dat$X)
    expect_identical(dat, original)
    expect_identical(.Random.seed, before)
})
