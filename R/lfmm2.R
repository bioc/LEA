# create class lfmm2
setClass("lfmm2Class",
         slots = c(K = "integer", 
                   lambda = "numeric",
                   B = "matrix",
                   U = "matrix",
                   V = "matrix"
                   )
)

lfmm2 <- function(input,
                  env, 
                  K, 
                  lambda = 1e-5,
                  effect.sizes = FALSE){

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
  if (is.character(env)){
    X <- read.env(env)
    if (anyNA(X)){
      stop("'env' file contains missing data (NA).")
    }
  } else {
    if (is.null(env)){
      stop("NULL value for argument 'env'.")
    }
    X <- as.matrix(env)
    if (anyNA(X)) {
      stop("The environmental matrix contains NA.")
    }
  }

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

# compute environmental effect sizes 
    if (effect.sizes){
    B <- (t(Ys - W) %*% Xs) %*% solve(t(Xs) %*% Xs + diag(lambda, nrow = d, ncol = d))
    B <- as.matrix(B)
    } else
    {B <-  matrix(NA)}
    
    

    obj <- new("lfmm2Class")
    obj@K <- as.integer(K)
    obj@lambda <- as.numeric(lambda)
    obj@B <- as.matrix(B)
    obj@U <- as.matrix(U)
    obj@V <- as.matrix(V)

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
            if (is.character(env)){
              X <- read.env(env)
              if (anyNA(X)){
                stop("'env' file contains missing data (NA).")
              }
            } else {
              if (is.null(env)){
                stop("NULL value for argument 'env'.")
              }
              X <- as.matrix(env)
              if (anyNA(X)) {
                stop("The environmental matrix contains NA.")
              }
            }
            
            d <-  ncol(X) #number of environmental variables
            n <-  nrow(X) #number of individuals
            
            if (nrow(Y) != n){
              stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
            }
            
            if (n < d) {
              stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
            }
   
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
              
              mod_lm =  lm(res_Y ~ res_X)
              statistics <- .lfmm2_full_statistics(mod_lm)
              r_squared <- statistics$adj.r.squared
              f_score <- statistics$fscores
              p_value <- statistics$pvalues
              
            } else {
            
            ## All p-values returned    
            if (linear){
              mod_lm <- lm(Y ~ ., data = data.frame(X, object@U)) 
              statistics <- .lfmm2_linear_statistics(mod_lm, d)
              p_value <- statistics$pvalues
              z_score <- statistics$zscores
            } else {
              na_action <- getOption("na.action")
              j <- 1L
              mod_glm <- glm(Y[, j] ~ ., data = data.frame(X, object@U), family = family)
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





