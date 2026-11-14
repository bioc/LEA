# create class lfmm2
setClass("lfmm2Class",
         slots = c(K = "integer", 
                   lambda = "numeric",
                   B = "matrix",
                   U = "matrix",
                   V = "matrix"
                   )
)

# Environmental preparation is deliberately separate from LFMM2 inference.
.lfmm2_numeric_frame <- function(env) {
    if (!is.data.frame(env) || !ncol(env) ||
        !all(vapply(env, function(x)
            (is.integer(x) || is.double(x)) && !is.object(x) &&
            is.null(dim(x)), logical(1))))
        stop("'env' must be a data.frame containing only ordinary numeric columns.")
}

# Validate language before evaluating any transformation. No user calls are run.
.lfmm2_formula_terms <- function(formula) {
    fail <- function() stop("Unsupported LFMM2 formula: use additive numeric variables, I(x^k) with a positive integer literal, log(x), or sqrt(x); interactions and arbitrary calls are not allowed.")
    if (!inherits(formula, "formula") || length(formula) != 2L) fail()
    rhs <- formula[[2L]]
    literal <- function(x, value) is.numeric(x) && length(x) == 1L &&
        is.finite(x) && x == value
    if (is.call(rhs) && identical(rhs[[1L]], as.name("-"))) {
        if (length(rhs) != 3L || !literal(rhs[[3L]], 1)) fail()
        rhs <- rhs[[2L]]
    }
    flatten <- function(x) {
        if (is.call(x) && identical(x[[1L]], as.name("+")) && length(x) == 3L)
            return(c(flatten(x[[2L]]), flatten(x[[3L]])))
        if (literal(x, 0)) return(list())
        list(x)
    }
    expressions <- flatten(rhs)
    if (!length(expressions)) fail()
    variable <- function(x) {
        if (!is.symbol(x) || identical(x, as.name("."))) fail()
        as.character(x)
    }
    terms <- lapply(expressions, function(x) {
        op <- "identity"; exponent <- NULL
        if (is.symbol(x)) {
            name <- variable(x)
        } else {
            if (!is.call(x) || length(x) != 2L || !is.symbol(x[[1L]])) fail()
            if (!is.null(names(x)) && any(nzchar(names(x)[-1L])))
                stop("Named arguments are not supported in LFMM2 formula transformations.")
            op <- as.character(x[[1L]])
            if (op == "I") {
                power <- x[[2L]]
                if (!is.call(power) || length(power) != 3L ||
                    !identical(power[[1L]], as.name("^"))) fail()
                if (!is.null(names(power)) && any(nzchar(names(power)[-1L])))
                    stop("Named arguments are not supported in LFMM2 formula transformations.")
                name <- variable(power[[2L]])
                exponent <- power[[3L]]
                if (!is.numeric(exponent) || length(exponent) != 1L ||
                    !is.finite(exponent) || exponent < 1 || exponent != floor(exponent)) fail()
                op <- "power"
            } else {
                if (!op %in% c("log", "sqrt")) fail()
                name <- variable(x[[2L]])
            }
        }
        list(variable = name, op = op, exponent = exponent,
             label = paste(deparse(x, width.cutoff = 500L, backtick = TRUE), collapse = ""))
    })
    labels <- vapply(terms, function(x) x$label, character(1))
    if (anyDuplicated(labels)) stop("Duplicated terms in LFMM2 formula.")
    list(terms = terms, labels = labels,
         normalized = paste("~", paste(labels, collapse = " + ")))
}

.lfmm2_env_spec <- function(object) {
    spec <- attr(object, "env.spec", exact = TRUE)
    if (is.null(spec)) return(NULL)
    fail <- function(detail)
        stop("Invalid or unsupported LFMM2 environmental metadata: ", detail, ".", call. = FALSE)
    if (typeof(spec) != "list" || is.object(spec) || !is.null(dim(spec)))
        fail("expected a named list")
    keys <- names(spec)
    required <- c("version", "formula", "mode", "predictors", "n", "d")
    if (is.null(keys) || anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys) ||
        length(keys) != length(required) || !all(required %in% keys))
        fail("missing, invalid or duplicated keys")
    if (!identical(spec[["version", exact = TRUE]], 1L)) fail("unsupported schema version")
    mode <- spec[["mode", exact = TRUE]]
    if (!is.character(mode) || !is.null(attributes(mode)) || length(mode) != 1L ||
        is.na(mode) || !mode %in% c("formula", "matrix", "data.frame", "file"))
        fail("invalid input mode")
    dimension <- function(x) is.integer(x) && is.null(attributes(x)) &&
        length(x) == 1L && !is.na(x) && x > 0L
    n <- spec[["n", exact = TRUE]]
    d <- spec[["d", exact = TRUE]]
    if (!dimension(n) || !dimension(d) || d > n) fail("invalid dimensions")
    predictors <- spec[["predictors", exact = TRUE]]
    if (!is.null(predictors) && (!is.character(predictors) ||
        !is.null(dim(predictors)) || length(predictors) != d))
        fail("invalid predictor structure")
    formula <- spec[["formula", exact = TRUE]]
    if (mode == "formula") {
        if (!is.character(formula) || !is.null(attributes(formula)) ||
            length(formula) != 1L || is.na(formula) || !nzchar(formula))
            fail("formula mode requires a formula string")
        if (d >= n || is.null(predictors) || anyNA(predictors) ||
            any(!nzchar(predictors)) || anyDuplicated(predictors))
            fail("formula dimensions or predictor structure")
        # Parse language only: malformed metadata must not execute R code.
        parsed <- tryCatch(parse(text = formula, keep.source = FALSE), error = function(e) NULL)
        if (length(parsed) != 1L || !is.call(parsed[[1L]]) ||
            !identical(parsed[[1L]][[1L]], as.name("~"))) fail("invalid formula")
        compiled <- tryCatch(.lfmm2_formula_terms(structure(parsed[[1L]],
            class = "formula", .Environment = baseenv())), error = function(e) NULL)
        if (is.null(compiled)) fail("invalid formula")
        if (length(compiled[["terms", exact = TRUE]]) != d ||
            !identical(compiled[["labels", exact = TRUE]], predictors))
            fail("formula dimensions or predictor structure")
    } else if (!is.null(formula)) {
        fail("formula must be NULL outside formula mode")
    }
    spec
}

.lfmm2_formula_diagnostic <- function(X) {
    # Subtract an anchor before centering to retain small variation around a
    # large absolute level. This copy never enters the statistical models.
    vapply(seq_len(ncol(X)), function(j) {
        x <- as.double(X[, j])
        shifted <- x - x[1L]
        if (any(!is.finite(shifted))) {
            x <- x / max(abs(x))
            shifted <- x - x[1L]
        }
        shifted <- shifted / max(abs(shifted))
        centered <- shifted - mean(shifted)
        centered / sqrt(sum(centered^2))
    }, numeric(nrow(X)))
}

.lfmm2_prepare_env <- function(env, formula = NULL, spec = NULL) {
    if (!is.null(spec) && !is.null(spec$formula)) {
        if (!is.character(spec$formula) || length(spec$formula) != 1L ||
            is.na(spec$formula)) stop("Invalid LFMM2 formula metadata.")
        formula <- stats::as.formula(spec$formula, env = baseenv())
    }
    if (is.null(formula)) {
        # Keep the historical file/matrix/vector path and its diagnostics.
        if (is.character(env)) {
            X <- read.env(env)
            if (anyNA(X)) stop("'env' file contains missing data (NA).")
            mode <- "file"
        } else {
            if (is.null(env)) stop("NULL value for argument 'env'.")
            if (is.data.frame(env)) .lfmm2_numeric_frame(env)
            X <- as.matrix(env)
            if (anyNA(X)) stop("The environmental matrix contains NA.")
            mode <- if (is.data.frame(env)) "data.frame" else "matrix"
        }
        normalized <- NULL
    } else {
        .lfmm2_numeric_frame(env)
        if (is.null(names(env)) || anyNA(names(env)) ||
            any(!nzchar(names(env))) || anyDuplicated(names(env)))
            stop("Formula 'env' columns must have unique, nonempty names.")
        compiled <- .lfmm2_formula_terms(formula)
        required <- vapply(compiled$terms, function(x) x$variable, character(1))
        missing <- setdiff(required, names(env))
        if (length(missing)) stop("Missing formula variables in 'env': ", paste(missing, collapse = ", "))
        if (!nrow(env)) stop("Formula 'env' must contain individuals.")
        values <- lapply(compiled$terms, function(term) {
            x <- env[[term$variable]]
            if (any(!is.finite(x))) stop("Non-finite values in formula variable '", term$variable, "'.")
            x <- suppressWarnings(switch(term$op, identity = x,
                power = x ^ term$exponent, log = log(x), sqrt = sqrt(x)))
            if (length(x) != nrow(env) || any(!is.finite(x)))
                stop("Formula term '", term$label, "' must produce finite numeric values for every individual.")
            if (all(x == x[1L])) stop("Constant formula predictor: ", term$label)
            x
        })
        X <- do.call(cbind, values)
        dimnames(X) <- list(rownames(env), compiled$labels)
        diagnostic <- .lfmm2_formula_diagnostic(X)
        if (qr(diagnostic, tol = .Machine$double.eps * max(dim(diagnostic)))$rank < ncol(X))
            stop("Linearly dependent formula predictors at numerical precision.")
        normalized <- compiled$normalized
        mode <- "formula"
        if (!is.null(spec) && (!identical(nrow(X), spec$n) ||
            !identical(ncol(X), spec$d) || !identical(colnames(X), spec$predictors)))
            stop("Formula environmental dimensions or predictor structure differ from the LFMM2 fit.")
    }
    list(X = X, spec = list(version = 1L, formula = normalized, mode = mode,
         predictors = colnames(X), n = nrow(X), d = ncol(X)))
}

# Check the fitted QR, not summary coefficient positions. For the uncommon
# rank-deficient case, put nuisance columns first to distinguish a redundant
# nuisance column from an environmental effect confounded with nuisance terms.
.lfmm2_check_estimable <- function(model, predictors, design = NULL,
                                  columns = seq_along(predictors) + 1L,
                                  expected = NULL) {
    retained <- model$qr$pivot[seq_len(model$rank)]
    bad <- !columns %in% retained
    if (!is.null(expected)) {
        labels <- if (is.matrix(model$coefficients)) rownames(model$coefficients)
                  else names(model$coefficients)
        bad <- bad | is.na(labels[columns]) | labels[columns] != expected
    }
    if (!any(bad) && model$rank < NROW(model$coefficients)) {
        if (is.null(design)) design <- stats::model.matrix(model)
        if (!is.null(model$weights)) {
            active <- model$weights > 0
            design <- design[active, , drop = FALSE] * sqrt(model$weights[active])
        }
        order <- c(setdiff(seq_len(ncol(design)), columns), columns)
        diagnostic <- qr(design[, order, drop = FALSE], tol = model$qr$tol)
        bad <- !columns %in% order[diagnostic$pivot[seq_len(diagnostic$rank)]]
    }
    if (any(bad))
        stop("Non-estimable environmental predictors in formula association model: ",
             paste(predictors[bad], collapse = ", "))
    invisible(NULL)
}

lfmm2 <- function(input,
                  env, 
                  K, 
                  lambda = 1e-5,
                  effect.sizes = FALSE,
                  formula = NULL){

    input.expression <- substitute(input)
    env.expression <- substitute(env)

## Check input response matrix 
## LEA  
    if (is.character(input)){
      Y <- read.lfmm(input)
      lst.unique <- unique(as.numeric(Y))
      if (9 %in% lst.unique){
        stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
      }
      if (-9 %in% lst.unique){
        stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
      }
    } else {
    ## Y is an R object       
      if (is.null(input)){
        stop("NULL value for argument 'input'.")
      }
      Y <- as.matrix(input)
      Y[Y == 9] <- NA
      Y[Y == -9] <- NA
      if (anyNA(Y)) {
        stop("The input matrix contains missing values: NA, 9 or -9 not allowed.")
      }
    }

## Check independent/covariate env matrix  
## LEA 
  prepared <- .lfmm2_prepare_env(env, formula = formula)
  X <- prepared$X

  if (length(K) > 1){
    stop("Multiple values of K not allowed.")
  }
  if (lambda <= 0){
    stop("The ridge regularization parameter must be positive.")
  }
  
  d <-  ncol(X) #number of environmental variables
  n <-  nrow(X) #number of individuals
  
  if (nrow(Y) != n){
    stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
  }
  
  if (n < d) {
    stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
  }
  
  fit.info <- list(version = 1L,
      genotype = list(label = .lfmm2_input_label(input.expression, input),
          mode = .lfmm2_input_mode(input), dimensions = dim(Y)),
      environment = list(label = .lfmm2_input_label(env.expression, env),
          mode = .lfmm2_input_mode(env),
          dimensions = if (is.character(env)) dim(X) else dim(as.matrix(env))),
      fitted.dimensions = dim(X), predictors = colnames(X),
      formula = prepared$spec$formula,
      predictor.sd = vapply(seq_len(ncol(X)), function(j) stats::sd(X[, j]), numeric(1)))
  names(fit.info$predictor.sd) <- colnames(X)

# centering  
  Xs <- scale(X, scale = FALSE)
  Ys <- scale(Y, scale = FALSE)

# run SVD of X: X = Q Sigma R
  
    svx <- svd(x = Xs, nu = n)
    Q <- svx$u
    
    d_lambda <- c(sqrt(lambda/(lambda + svx$d^2)), rep(1, n-d))
    d_lambda_inv <- c(sqrt((lambda + svx$d^2)/lambda), rep(1, n-d))
    D_inv <- diag(d_lambda_inv)
    D  <- diag(d_lambda)
 
# run SVD of modified Y    
    svk <- svd(D %*% t(Q) %*% Ys, nu = K)

    if (K > 1) {
      Sigma_k <- diag(svk$d[1:K])
    } else {
      Sigma_k <- as.matrix(svk$d[1])
    }
    
# compute the latent matrix W
    W <- Q %*% D_inv %*% tcrossprod(svk$u %*% Sigma_k, svk$v[,1:K])

# compute LFMM factors U and loadings V
# Non orthogonal factors
    U <- crossprod(t(Q %*% D_inv), svk$u %*% Sigma_k)
    #U <- Q %*% D_inv %*% svk$u %*% Sigma_k
    V <- svk$v[,1:K]
    # The full SVD output is no longer needed during coefficient computation.
    rm(svk)

# compute environmental effect sizes 
    if (effect.sizes){
    B <- (t(Ys - W) %*% Xs) %*% solve(t(Xs) %*% Xs + diag(lambda, nrow = d, ncol = d))
    B <- as.matrix(B)
    Cb <- crossprod(B) / nrow(B)
    } else {
        # Use the identical ridge coefficient definition in bounded locus blocks.
        ridge <- solve(t(Xs) %*% Xs + diag(lambda, nrow = d, ncol = d))
        Cb <- matrix(0, d, d)
        block.size <- max(1L, floor(1e6 / max(n, d)))
        for (first in seq.int(1L, ncol(Ys), by = block.size)) {
            loci <- seq.int(first, min(ncol(Ys), first + block.size - 1L))
            coefficients <- crossprod(Ys[, loci, drop = FALSE] - W[, loci, drop = FALSE], Xs) %*% ridge
            Cb <- Cb + crossprod(coefficients)
        }
        Cb <- Cb / ncol(Ys)
        B <- matrix(NA)
    }
    dimnames(Cb) <- list(colnames(X), colnames(X))
    
    

    obj <- new("lfmm2Class")
    obj@K <- as.integer(K)
    obj@lambda <- as.numeric(lambda)
    obj@B <- as.matrix(B)
    obj@U <- as.matrix(U)
    obj@V <- as.matrix(V)
    attr(obj, "env.spec") <- prepared$spec
    attr(obj, "fit.info") <- fit.info
    attr(obj, "gea.metric") <- list(version = 1L, Cb = Cb)

## LEA 
    return(obj)
}


            
# Fast summaries are restricted to regular multivariate fits. Keep summary.lm's
# historical errors, alias extraction and warnings for exceptional fits.
.lfmm2_regular_fit <- function(model) {
    if (!inherits(model, "mlm") || model$rank != nrow(model$coefficients) ||
        model$df.residual <= 0 || model$rank <= 1L ||
        anyNA(colnames(model$coefficients))) return(FALSE)
    variance <- colSums(model$residuals^2) / model$df.residual
    if (any(!is.finite(variance))) return(FALSE)
    # Conservative upper bound for mean(f)^2 + var(f): delegate every possible
    # "essentially perfect fit" warning to the original summaries.
    bound <- colSums(model$fitted.values^2) /
        (nrow(model$fitted.values) - 1L)
    !any(variance < bound * 1e-30)
}

.lfmm2_response_names <- function(model) {
    labels <- colnames(model$coefficients)
    defaults <- paste0("Y", seq_len(ncol(model$coefficients)))
    if (is.null(labels)) labels <- defaults
    labels[labels == ""] <- defaults[labels == ""]
    paste("Response", labels)
}

.lfmm2_linear_statistics <- function(model, d) {
    if (!.lfmm2_regular_fit(model) || d < 1L) {
        sm <- summary(model)
        return(list(pvalues = sapply(sm, function(x) x$coeff[2:(d + 1), 4]),
                    zscores = as.matrix(sapply(sm, function(x) x$coeff[2:(d + 1), 3]))))
    }
    rank <- model$rank
    selected <- seq_len(d) + 1L
    # summary.lm uses the pivot order, including its historical positional
    # coefficient extraction. Rank-deficient fits use the fallback above.
    coefficients <- model$coefficients[model$qr$pivot[seq_len(rank)], , drop = FALSE]
    unscaled <- diag(chol2inv(model$qr$qr[seq_len(rank), seq_len(rank), drop = FALSE]))
    variance <- colSums(model$residuals^2) / model$df.residual
    statistic <- coefficients[selected, , drop = FALSE] /
        sqrt(outer(unscaled[selected], variance))
    probability <- 2 * stats::pt(abs(statistic), model$df.residual, lower.tail = FALSE)
    labels <- .lfmm2_response_names(model)
    if (d == 1L) {
        probability <- as.vector(probability)
        names(probability) <- labels
        statistic <- matrix(as.vector(statistic), ncol = 1L,
                            dimnames = list(labels, NULL))
    } else {
        dimnames(statistic) <- dimnames(probability) <-
            list(rownames(coefficients)[selected], labels)
    }
    list(pvalues = probability, zscores = statistic)
}

.lfmm2_full_statistics <- function(model) {
    if (!.lfmm2_regular_fit(model)) {
        sm <- summary(model)
        return(list(adj.r.squared = sapply(sm, function(x) x$adj.r.squared),
                    fscores = sapply(sm, function(x) x$fstat[1]),
                    pvalues = sapply(sm, function(x)
                        pf(x$fstat[1], x$fstat[2], x$fstat[3], lower.tail = FALSE))))
    }
    fitted <- model$fitted.values
    mss <- colSums(sweep(fitted, 2L, colMeans(fitted))^2)
    rss <- colSums(model$residuals^2)
    rdf <- model$df.residual
    # Preserve summary.lm's historical second-stage degrees of freedom, NOT
    # those of a conventional partial F test including the nuisance factors.
    numerator_df <- model$rank - 1L
    statistic <- (mss / numerator_df) / (rss / rdf)
    adjusted <- 1 - (1 - mss / (mss + rss)) *
        ((nrow(fitted) - 1L) / rdf)
    labels <- .lfmm2_response_names(model)
    names(statistic) <- paste0(labels, ".value")
    names(adjusted) <- labels
    list(adj.r.squared = adjusted, fscores = statistic,
         pvalues = pf(statistic, numerator_df, rdf, lower.tail = FALSE))
}

# Custom family factories, initializers and callbacks retain formula-based glm.
.lfmm2_glm_reusable_family <- function(family) {
    if (!is.list(family) || !inherits(family, "family") ||
        !is.character(family$family) || length(family$family) != 1L ||
        anyNA(family$family)) return(FALSE)
    factory <- switch(family$family,
        gaussian = stats::gaussian, binomial = stats::binomial,
        poisson = stats::poisson, Gamma = stats::Gamma,
        inverse.gaussian = stats::inverse.gaussian,
        quasibinomial = stats::quasibinomial, quasipoisson = stats::quasipoisson)
    if (is.null(factory) || !is.character(family$link) || length(family$link) != 1L)
        return(FALSE)
    reference <- tryCatch(factory(link = family$link), error = function(e) NULL)
    if (is.null(reference) || !identical(family$initialize, reference$initialize))
        return(FALSE)
    callbacks <- family[vapply(family, is.function, logical(1))]
    all(vapply(callbacks, function(fun)
        !is.null(environment(fun)) &&
        identical(topenv(environment(fun)), asNamespace("stats")), logical(1)))
}

.lfmm2_glm_reusable_na_action <- function(action, envir) {
    if (is.null(action)) return(TRUE)
    if (is.character(action) && length(action) == 1L && !anyNA(action)) {
        action <- tryCatch(get(action, envir = envir, mode = "function"),
                           error = function(e) NULL)
    }
    if (!is.function(action)) return(FALSE)
    standard <- list(stats::na.omit, stats::na.exclude, stats::na.fail, stats::na.pass)
    any(vapply(standard, function(fun) identical(action, fun), logical(1)))
}

.lfmm2_glm_refit <- function(design, response, family, template) {
    frame_response <- response
    names(response) <- names(template$y)
    model <- stats::glm.fit(x = design, y = response,
        weights = NULL, start = NULL, etastart = NULL, mustart = NULL,
        offset = NULL, family = family, control = template$control,
        intercept = attr(template$terms, "intercept") > 0L, singular.ok = TRUE)
    # Preserve glm's summary metadata and raw model response, while sharing the
    # fixed predictor columns. IRLS, dispersion and inference still use stats.
    model$model <- template$model
    names(frame_response) <- names(template$model[[1L]])
    model$model[[1L]] <- frame_response
    fields <- c("call", "formula", "terms", "data", "offset", "control",
                "method", "contrasts", "xlevels")
    model[fields] <- template[fields]
    class(model) <- c("glm", "lm")
    model
}

setGeneric("lfmm2.test", function(object, input, env, 
                                  full = FALSE,
                                  genomic.control = TRUE, 
                                  linear = TRUE, 
                                  family = binomial(link = "logit")) matrix);
setMethod("lfmm2.test", "lfmm2Class",
          function(object, 
                   input,
                   env,
                   full,
                   genomic.control, 
                   linear,
                   family
                   ) {

            ## Check input matrix   
            ## LEA  
            if (is.character(input)){
              warning("Reading large input files with 'read.lfmm()' may be slow. See 'data.table::fread()' for fast import.")
              Y <- read.lfmm(input)
              lst.unique <- unique(as.numeric(Y))
              if (9 %in% lst.unique){
                stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
              }
              if (-9 %in% lst.unique){
                stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
              }
            } else {
              ## Y is an R object       
              if (is.null(input)){
                stop("NULL value for argument 'input'.")
              }
              Y <- as.matrix(input)
              Y[Y == 9] <- NA
              Y[Y == -9] <- NA
              if (anyNA(Y)) {
                stop("The input matrix contains missing values (NA or 9).")
              }
            }
            
            ## Check independent/covariate matrix  
            ## LEA 
            spec <- .lfmm2_env_spec(object)
            X <- .lfmm2_prepare_env(env, spec = spec)$X
            formula_fit <- !is.null(spec) && !is.null(spec$formula)
            if (formula_fit && nrow(object@U) != nrow(X))
                stop("Formula individuals do not match the LFMM2 latent factors.")

            d <-  ncol(X) #number of environmental variables
            n <-  nrow(X) #number of individuals
            
            if (nrow(Y) != n){
              stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
            }
            
            if (n < d) {
              stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
            }
   
            environment_names <- if (formula_fit) names(data.frame(X, object@U))[seq_len(d)] else NULL
            p <- ncol(Y)
            gif <-  NULL
            p_value <- NULL
            z_score <- NULL
            f_score <- NULL           
            r_squared <- NULL
            
            if (full){
              
              ## a single p-value is returned (f-test)
              ## Check linear models  
              if (linear == FALSE){
                stop("Option full == TRUE is available only for linear models.")
              }
              
              ## partial regression 
              mod_Y <- lm(Y ~ ., data = data.frame(object@U)) 
              res_Y = mod_Y$residuals
              if (p > 1L && is.null(mod_Y$na.action) &&
                  mod_Y$rank == ncol(object@U) + 1L) {
                # Reuse the same nuisance QR for X. With omitted rows or an
                # aliased nuisance design, retain the original lm() path.
                res_X <- qr.resid(mod_Y$qr, X)
                if (d == 1L) res_X <- drop(res_X)
              } else {
                mod_X = lm(X ~ ., data = data.frame(object@U))
                res_X = mod_X$residuals
              }
              
              if (formula_fit) {
                # Use the nuisance columns retained by the actual first stage.
                # A normalized diagnostic detects absorbed effects without
                # rejecting a valid residual regression due to raw X levels.
                nuisance <- stats::model.matrix(mod_Y)
                rows <- match(rownames(mod_Y$model), rownames(data.frame(object@U)))
                if (anyNA(rows)) stop("Cannot align formula predictors with the global association model.")
                nuisance <- nuisance[, mod_Y$qr$pivot[seq_len(mod_Y$rank)], drop = FALSE]
                joint <- cbind(nuisance, .lfmm2_formula_diagnostic(X)[rows, , drop = FALSE])
                joint_fit <- stats::lm.fit(joint, numeric(nrow(joint)),
                    tol = .Machine$double.eps * max(dim(joint)))
                .lfmm2_check_estimable(joint_fit, spec$predictors, design = joint,
                    columns = ncol(nuisance) + seq_len(d))
              }
              mod_lm =  lm(res_Y ~ res_X)
              if (formula_fit) .lfmm2_check_estimable(mod_lm, spec$predictors)
              statistics <- .lfmm2_full_statistics(mod_lm)
              r_squared <- statistics$adj.r.squared
              f_score <- statistics$fscores
              p_value <- statistics$pvalues
              
            } else {
            
            ## All p-values returned    
            if (linear){
              mod_lm <- lm(Y ~ ., data = data.frame(X, object@U))
              if (formula_fit) .lfmm2_check_estimable(mod_lm, spec$predictors,
                                                     expected = environment_names)
              statistics <- .lfmm2_linear_statistics(mod_lm, d)
              p_value <- statistics$pvalues
              z_score <- statistics$zscores
            } else {
              na_action <- getOption("na.action")
              j <- 1L
              mod_glm <- glm(Y[, j] ~ ., data = data.frame(X, object@U), family = family)
              if (formula_fit) .lfmm2_check_estimable(mod_glm, spec$predictors,
                                                     expected = environment_names)
              sm <- summary(mod_glm)
              probability <- sm$coeff[2:(d + 1), 4]
              statistic <- sm$coeff[2:(d + 1), 3]
              # Match rbind()'s first-row names and singleton simplification.
              # Allocate once, before iterating over the remaining loci.
              p_value <- matrix(NA_real_, p, length(probability))
              z_score <- matrix(NA_real_, p, length(statistic))
              if (!is.null(names(probability))) colnames(p_value) <- names(probability)
              if (!is.null(names(statistic))) colnames(z_score) <- names(statistic)
              p_value[1L, ] <- probability
              z_score[1L, ] <- statistic
              # The first ordinary glm supplies exact terms, contrasts,
              # observation names and controls for subsequent fits.
              reuse_design <- p > 1L &&
                  "singular.ok" %in% names(formals(stats::glm.fit)) &&
                  is.numeric(Y) && is.numeric(X) &&
                  is.numeric(object@U) && is.null(attr(Y, "class")) &&
                  is.null(attr(X, "class")) && is.null(attr(object@U, "class")) &&
                  all(is.finite(Y)) && is.null(mod_glm$na.action) &&
                  nrow(mod_glm$model) == n &&
                  identical(na_action, getOption("na.action")) &&
                  .lfmm2_glm_reusable_na_action(na_action, environment(mod_glm$formula)) &&
                  is.null(mod_glm$offset) && .lfmm2_glm_reusable_family(family)
              if (reuse_design) {
                template <- mod_glm
                design <- stats::model.matrix(template$terms, template$model,
                                              contrasts.arg = template$contrasts)
              }
              for (j in seq_len(p)[-1L]) {
                if (reuse_design) {
                  mod_glm <- .lfmm2_glm_refit(design, Y[, j], family, template)
                } else {
                  mod_glm <- glm(Y[, j] ~ ., data = data.frame(X, object@U), family = family)
                }
                if (formula_fit) .lfmm2_check_estimable(mod_glm, spec$predictors,
                                                       expected = environment_names)
                sm <- summary(mod_glm)
                p_value[j, ] <- sm$coeff[2:(d + 1), 4]
                z_score[j, ] <- sm$coeff[2:(d + 1), 3]
              }
            }
            
            }
              
              
            if (genomic.control){
              
            if (!full){
              
              if (d == 1){
                gif <- median(z_score^2)/qchisq(0.5, df = 1, lower.tail = FALSE)
              } else if (linear) {
                gif <- apply(z_score^2, 1, median)/qchisq(0.5, df = 1, lower.tail = FALSE)
              } else {
                # GLM statistics have loci in rows and covariates in columns.
                gif <- apply(z_score^2, 2, median)/qchisq(0.5, df = 1, lower.tail = FALSE)
              }
              if (!linear && d > 1) {
                p_value <- pchisq(sweep(z_score^2, 2, gif, "/"), df = 1, lower.tail = FALSE)
              } else {
                p_value <- pchisq(z_score^2/gif, df = 1, lower.tail = FALSE)
              }
                
            } else {
              
              gif <- median(f_score)/qf(0.5, d, n - d - 1, lower.tail = FALSE)
              p_value <- pf(f_score/gif, d, n - d - 1, lower.tail = FALSE)
            }
            
            }
            
            if (!full & anyNA(z_score)) {
              warning("NA in significance values and z-scores. Check the input matrix for non-variable genomic sites.")
            }
            if (full & anyNA(r_squared)) {
              warning("NA in significance values and R-squared. Check the input matrix for non-variable genomic sites.")
            }
            
       res <- list(pvalues = p_value, zscores = z_score, fscores = f_score, adj.r.squared = r_squared,  gif = gif)
       return(res)
      }
)





