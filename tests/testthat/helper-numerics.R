# Deterministic fixtures restore the caller's RNG state and use explicit RNGs.
make_lfmm_fixture <- function(n, p, d, K, seed, binary = FALSE) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    old_kind <- RNGkind()
    on.exit({
        do.call(RNGkind, as.list(old_kind))
        if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
    })
    # No sample() calls: sample.kind is irrelevant and was unavailable in R 3.3.
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
    X <- matrix(rnorm(n * d), n, d)
    A <- matrix(seq(0.15, 0.45, length.out = d * K), d, K)
    U <- X %*% A + matrix(rnorm(n * K), n, K)
    V <- sweep(matrix(rnorm(p * K), p, K), 2, seq(1.5, 0.7, length.out = K), "*")
    B <- matrix(0, p, d)
    B[seq_len(min(6L, p)), ] <- matrix(
        seq(-0.8, 0.9, length.out = min(6L, p) * d), min(6L, p), d)
    signal <- X %*% t(B) + tcrossprod(U, V)
    if (binary) {
        prob <- plogis(0.25 * signal)
        Y <- matrix(rbinom(n * p, 1, prob), n, p)
    } else Y <- signal + matrix(rnorm(n * p, sd = 0.7), n, p)
    list(X = X, Y = Y, K = K, B = B)
}

lfmm_fixtures <- list(
    simple = make_lfmm_fixture(60L, 40L, 1L, 1L, 7101L),
    multi = make_lfmm_fixture(80L, 120L, 3L, 2L, 7102L),
    binary = make_lfmm_fixture(80L, 40L, 3L, 2L, 7103L, binary = TRUE)
)

# Numeric values and missing/nonfinite patterns are checked separately from
# attributes. Explicit elementwise absolute + relative bounds avoid a global
# scale hiding errors in small entries.
expect_numeric_close <- function(actual, expected, rtol = 1e-8, atol = 1e-9) {
    a <- as.numeric(actual)
    b <- as.numeric(expected)
    expect_length(a, length(b))
    expect_identical(is.na(a), is.na(b))
    expect_identical(is.nan(a), is.nan(b))
    expect_identical(is.infinite(a), is.infinite(b))
    if (length(a) != length(b)) return(invisible(NULL))
    finite <- is.finite(a) & is.finite(b)
    expect_true(all(abs(a[finite] - b[finite]) <= atol + rtol * abs(b[finite])),
                info = "elementwise absolute + relative numerical tolerance")
    inf <- is.infinite(b)
    expect_identical(a[inf], b[inf])
}

expect_p_close <- function(actual, expected) {
    expect_numeric_close(actual, expected, rtol = 1e-7, atol = 1e-12)
    a <- as.numeric(actual)
    b <- as.numeric(expected)
    small <- is.finite(b) & b > 0 & b < 1e-8
    if (any(small)) {
        expect_true(all(a[small] > 0))
        expect_numeric_close(log(a[small]), log(b[small]),
                             rtol = 1e-7, atol = 1e-9)
    }
    zero <- is.finite(b) & b == 0
    expect_identical(a[zero], b[zero])
}

expect_result_close <- function(actual, expected, glm = FALSE) {
    expect_identical(names(actual), names(expected))
    for (nm in names(expected)) {
        expect_identical(is.null(actual[[nm]]), is.null(expected[[nm]]))
        expect_identical(attributes(actual[[nm]]), attributes(expected[[nm]]))
        if (nm == "pvalues") expect_p_close(actual[[nm]], expected[[nm]])
        else expect_numeric_close(actual[[nm]], expected[[nm]],
                                  rtol = if (glm) 1e-6 else 1e-8,
                                  atol = if (glm) 1e-7 else 1e-9)
    }
}

historical_test <- function(fit, Y, X, full = FALSE, genomic.control = TRUE,
                            linear = TRUE, family = binomial()) {
    historical_lfmm2_test(fit, Y, X, full, genomic.control, linear, family)
}

# Independent symmetric ridge-filter reference, not the production SVD basis.
# Restrict to moderate lambda: I - ridge hat matrix suffers cancellation when
# lambda is tiny. The frozen oracle covers the default 1e-5 separately.
reference_fit <- function(Y, X, K, lambda) {
    Xc <- scale(X, scale = FALSE)
    Yc <- scale(Y, scale = FALSE)
    gram <- crossprod(Xc) + diag(lambda, ncol(Xc))
    H <- diag(nrow(Xc)) - Xc %*% solve(gram, t(Xc))
    eig <- eigen(H, symmetric = TRUE)
    stopifnot(all(eig$values > 0))
    root <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
    invroot <- eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)
    s <- svd(root %*% Yc)
    keep <- seq_len(K)
    W <- invroot %*% sweep(s$u[, keep, drop = FALSE], 2, s$d[keep], "*") %*%
        t(s$v[, keep, drop = FALSE])
    B <- t(solve(gram, crossprod(Xc, Yc - W)))
    list(W = W, B = B)
}

# Per-locus lm/glm references avoid the production multivariate-summary path.
reference_association <- function(Y, X, U, family = NULL) {
    d <- ncol(X)
    design <- data.frame(X, U)
    models <- lapply(seq_len(ncol(Y)), function(j) {
        y <- Y[, j]
        if (is.null(family)) lm(y ~ ., data = design)
        else glm(y ~ ., data = design, family = family)
    })
    coef <- lapply(models, function(m) summary(m)$coefficients)
    list(stat = vapply(coef, function(c) c[seq_len(d) + 1L, 3], numeric(d)),
         p = vapply(coef, function(c) c[seq_len(d) + 1L, 4], numeric(d)),
         models = models)
}

# Historical two-stage full test, computed per locus.
reference_full_historical <- function(Y, X, U) {
    nuisance <- cbind(1, U)
    ry <- qr.resid(qr(nuisance), Y)
    rx <- qr.resid(qr(nuisance), X)
    summaries <- lapply(seq_len(ncol(Y)), function(j) summary(lm(ry[, j] ~ rx)))
    f <- vapply(summaries, function(s) unname(s$fstatistic[1]), numeric(1))
    p <- vapply(summaries, function(s) {
        fs <- s$fstatistic
        pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
    }, numeric(1))
    adj <- vapply(summaries, function(s) s$adj.r.squared, numeric(1))
    list(f = f, p = p, adj = adj)
}

# Conventional nested-model partial F reference uses effective model ranks.
# It does not claim any improvement in genomic-control calibration or power.
reference_partial_f <- function(Y, X, U) {
    reduced <- cbind(1, U)
    complete <- cbind(reduced, X)
    q0 <- qr(reduced)
    q1 <- qr(complete)
    df1 <- q1$rank - q0$rank
    df2 <- nrow(Y) - q1$rank
    stopifnot(df1 > 0, df2 > 0)
    rss0 <- colSums(qr.resid(q0, Y)^2)
    rss1 <- colSums(qr.resid(q1, Y)^2)
    f <- ((rss0 - rss1) / df1) / (rss1 / df2)
    list(f = f, p = pf(f, df1, df2, lower.tail = FALSE),
         df1 = df1, df2 = df2)
}
