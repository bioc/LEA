# Optional attributes leave the historical five-slot serialization intact.
.lfmm2_input_mode <- function(x) {
    if (is.character(x)) "file" else if (is.data.frame(x)) "data.frame" else "matrix"
}
.lfmm2_input_label <- function(expression, value) {
    if (is.character(value) && length(value) == 1L) return(value)
    label <- paste(deparse(expression, width.cutoff = 60L), collapse = " ")
    if (nchar(label) > 160L) label <- paste0(substr(label, 1L, 157L), "...")
    label
}
.lfmm2_retained <- function(object) {
    B <- object@B
    is.numeric(B) && length(B) > 0L && all(is.finite(B)) &&
        nrow(B) == nrow(object@V)
}
.lfmm2_fit_info <- function(object) {
    info <- attr(object, "fit.info", exact = TRUE)
    if (!is.list(info) || !identical(info$version, 1L)) return(NULL)
    dimensions <- function(x) is.numeric(x) && length(x) == 2L &&
        all(is.finite(x)) && all(x > 0) && all(x == floor(x))
    input <- function(x) is.list(x) && is.character(x$label) &&
        length(x$label) == 1L && !is.na(x$label) && dimensions(x$dimensions)
    if (!input(info$genotype) || !input(info$environment) ||
        !dimensions(info$fitted.dimensions)) return(NULL)
    if (!is.null(info$formula) && (!is.character(info$formula) ||
        length(info$formula) != 1L || is.na(info$formula))) return(NULL)
    info
}
.lfmm2_description <- function(object) {
    info <- .lfmm2_fit_info(object)
    spec <- tryCatch(.lfmm2_env_spec(object), error = function(e) NULL)
    describe <- function(x) {
        if (is.null(x)) return("unavailable (provenance not recorded)")
        paste0(x$label, " (", paste(x$dimensions, collapse = " x "), ")")
    }
    list(genotypes = describe(info$genotype), environment = describe(info$environment),
         fit.info = info, K = object@K, lambda = object@lambda,
         predictors = if (!is.null(info)) info$fitted.dimensions[2L] else
             if (.lfmm2_retained(object)) ncol(object@B) else
             if (!is.null(spec)) spec$d else NA_integer_,
         coefficients.retained = .lfmm2_retained(object),
         coefficient.dimensions = if (.lfmm2_retained(object)) dim(object@B) else NULL,
         formula = if (!is.null(info)) info$formula else spec$formula, latent.factors = "object@U")
}
.lfmm2_show <- function(object) {
    x <- .lfmm2_description(object)
    cat("LFMM2 fit\nGenotypes: ", x$genotypes, "\nEnvironment: ", x$environment,
        "\nLatent factors: K = ", x$K, " (object@U)\nGEA coefficients: ",
        if (x$coefficients.retained) paste0("retained (", paste(x$coefficient.dimensions, collapse = " x "), ")") else "not retained",
        "\nUse summary(object) for GEA spectral diagnostics.\n", sep = "")
    invisible(object)
}
setMethod("show", "lfmm2Class", function(object) .lfmm2_show(object))
setMethod("print", "lfmm2Class", function(x, ...) .lfmm2_show(x))

.lfmm2_spectrum <- function(object) {
    unavailable <- function(reason) list(available = FALSE, reason = reason)
    metric <- attr(object, "gea.metric", exact = TRUE)
    source <- "gea.metric"
    if (is.null(metric)) {
        if (!.lfmm2_retained(object)) return(unavailable("neither Cb nor valid retained coefficients are available"))
        Cb <- crossprod(object@B) / nrow(object@B)
        source <- "retained B"
    } else {
        if (!is.list(metric) || !identical(metric$version, 1L))
            return(unavailable("invalid or unsupported gea.metric metadata"))
        Cb <- metric$Cb
    }
    if (!is.matrix(Cb) || !is.numeric(Cb) || !nrow(Cb) ||
        nrow(Cb) != ncol(Cb) || any(!is.finite(Cb)))
        return(unavailable("Cb must be a finite numeric square matrix"))
    d <- nrow(Cb)
    info <- .lfmm2_fit_info(object)
    if (!is.null(info) && d != info$fitted.dimensions[2L])
        return(unavailable("Cb dimensions disagree with fitted predictors"))
    scale <- max(abs(Cb))
    # Work at unit scale to avoid overflow and retain tiny nonzero metrics.
    normalized <- if (scale > 0) Cb / scale else Cb
    relative.tolerance <- 100 * .Machine$double.eps * d
    tolerance <- relative.tolerance * scale
    if (max(abs(normalized - t(normalized))) > relative.tolerance)
        return(unavailable("Cb is not symmetric within numerical tolerance"))
    eig <- eigen((normalized + t(normalized)) / 2, symmetric = TRUE)
    if (min(eig$values) < -relative.tolerance)
        return(unavailable("Cb has materially negative eigenvalues"))
    values <- pmax(eig$values, 0)
    total <- sum(values)
    proportion <- cumulative <- rep(NA_real_, d)
    axes99 <- NA_integer_
    if (total > 0) {
        proportion <- values / total
        cumulative <- cumsum(proportion)
        axes99 <- which(cumulative >= 0.99 - relative.tolerance)[1L]
    }
    if (any(!is.finite(values * scale)))
        return(unavailable("Cb eigenvalues exceed the finite numeric range"))
    labels <- colnames(Cb)
    if (is.null(labels) && !is.null(info) && is.character(info$predictors) &&
        length(info$predictors) == d) labels <- info$predictors
    if (is.null(labels)) labels <- paste0("Predictor", seq_len(d))
    rownames(eig$vectors) <- labels
    colnames(eig$vectors) <- paste0("Axis", seq_len(d))
    list(available = TRUE, source = source, eigenvalues = values * scale,
         vectors = eig$vectors, proportion = proportion, cumulative = cumulative,
         axes99 = axes99, variables = labels, tolerance = tolerance,
         relative.tolerance = relative.tolerance)
}
setMethod("summary", "lfmm2Class", function(object, ...) {
    result <- .lfmm2_description(object)
    result$gea.spectrum <- .lfmm2_spectrum(object)
    class(result) <- "summary.lfmm2"
    result
})
print.summary.lfmm2 <- function(x, ...) {
    cat("LFMM2 model summary\n\nGenotypes:      ", x$genotypes,
        "\nEnvironment:    ", x$environment,
        "\nFitted predictors: ", x$predictors,
        "\nLatent factors: K = ", x$K, " (object@U)",
        "\nRegularization: lambda = ", format(x$lambda), "\n", sep = "")
    if (!is.null(x$formula)) cat("Formula: ", x$formula, "\n", sep = "")
    cat("\nGEA coefficients:\n  Per-locus coefficients ",
        if (x$coefficients.retained) "retained" else "not retained",
        "\n\nGEA coefficient spectrum (Cb):\n", sep = "")
    spectrum <- x$gea.spectrum
    if (!spectrum$available) cat("  Unavailable: ", spectrum$reason, ".\n", sep = "")
    else if (is.na(spectrum$axes99)) cat("  Zero trace: spectral proportions and the 99% axis count are undefined.\n")
    else {
        d <- length(spectrum$eigenvalues)
        keep <- seq_len(if (d > 12L) spectrum$axes99 else d)
        table <- rbind("Contribution (%)" = 100 * spectrum$proportion[keep],
                       "Cumulative (%)" = 100 * spectrum$cumulative[keep])
        colnames(table) <- paste("Axis", keep)
        print(round(table, 1), quote = FALSE)
        if (length(keep) < d) cat(d - length(keep), "additional axes omitted; full spectrum is available in the summary result.\n")
        cat("\n", spectrum$axes99, " axes account for 99% of the GEA coefficient spectrum.\n", sep = "")
    }
    cat("\nNote: Standardizing environmental predictors to unit\nstandard deviation (sd = 1) is recommended before fitting\nLFMM2 for meaningful comparisons of spectral contributions.\n")
    invisible(x)
}
