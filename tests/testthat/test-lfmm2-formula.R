formula_env <- function() {
    dat <- lfmm_fixtures$binary
    data.frame(x = dat$X[, 1], y = dat$X[, 2],
               z = dat$X[, 3], w = exp(dat$X[, 3]/4) + 1)
}

formula_manual <- function(env, formula) {
    # Independent reference for the supported numeric subset, not a production helper.
    mm <- model.matrix(formula, model.frame(formula, env, na.action = na.pass))
    mm[, attr(mm, "assign") != 0L, drop = FALSE]
}

test_that("formula preparation matches manual transformations without centering or intercept", {
    env <- formula_env()
    formulas <- list(~ x + y + z, ~ x + I(x^2), ~ x + I(x^3),
                     ~ log(w), ~ sqrt(w), ~ 0 + y + x, ~ y + x - 1,
                     ~ x + y + I(z^2) + log(w))
    for (f in formulas) {
        actual <- LEA:::.lfmm2_prepare_env(env, f)$X
        manual <- formula_manual(env, f)
        expect_identical(actual, manual)
        expect_false("(Intercept)" %in% colnames(actual))
    }
    power <- LEA:::.lfmm2_prepare_env(data.frame(x = seq_len(10) + 10), ~ I(x^2))$X
    expect_identical(as.vector(power), (seq_len(10) + 10)^2)
    expect_equal(as.vector(scale(power, scale = FALSE)),
                 (seq_len(10) + 10)^2 - mean((seq_len(10) + 10)^2))
    expect_false(isTRUE(all.equal(as.vector(scale(power, scale=FALSE)),
                                 as.vector(scale(seq_len(10) + 10, scale=FALSE))^2)))
})

test_that("formula fit and association preserve manual matrix inference in every mode", {
    env <- formula_env(); Y <- lfmm_fixtures$binary$Y
    for (f in list(~ x, ~ y + x + I(z^2) + log(w), ~ x + I(x^3) + sqrt(w))) {
        X <- formula_manual(env, f)
        for (effects in c(FALSE, TRUE)) {
            fit <- lfmm2(Y, env, 2L, 0.1, effects, formula=f)
            manual <- lfmm2(Y, X, 2L, 0.1, effects)
            for (slot in c("K", "lambda", "B", "U", "V"))
                expect_identical(slot(fit, slot), slot(manual, slot))
            for (mode in c("linear", "full", "glm")) for (control in c(FALSE, TRUE)) {
                actual <- lfmm2.test(fit, Y, env, full=mode=="full",
                                     linear=mode!="glm", genomic.control=control)
                expected <- lfmm2.test(manual, Y, X, full=mode=="full",
                                       linear=mode!="glm", genomic.control=control)
                expect_identical(actual, expected)
            }
        }
    }
    expect_false("formula" %in% names(formals(lfmm2.test)))
    expect_identical(tail(names(formals(lfmm2)), 1), "formula")
})

test_that("formula reconstruction uses stored names and order, not session bindings", {
    env <- formula_env(); Y <- lfmm_fixtures$binary$Y
    f <- local({ hidden <- new.env(); hidden$large <- rep(0,10000)
                 structure(~ y + I(x^2), .Environment=hidden) })
    fit <- lfmm2(Y, env, 2L, formula=f)
    spec <- attr(fit, "env.spec", exact=TRUE)
    expect_type(spec$formula, "character")
    expect_null(attr(spec$formula, ".Environment"))
    expect_identical(spec$predictors, c("y", "I(x^2)"))
    expected <- lfmm2.test(fit,Y,env)
    expect_identical(lfmm2.test(fit,Y,env[,c("w","x","z","y")]), expected)
    env$unused <- 1
    expect_identical(lfmm2.test(fit,Y,env), expected)
    expect_error(lfmm2.test(fit,Y,env[,setdiff(names(env),"x")]), "Missing formula variables")
    expect_error(lfmm2.test(fit,Y[-1,],env[-1,]), "dimensions or predictor structure")
    expect_error(lfmm2.test(fit,Y,as.matrix(env)), "data.frame")
    broken <- fit; s <- spec; s$predictors <- rev(s$predictors); attr(broken,"env.spec") <- s
    expect_error(lfmm2.test(broken,Y,env), "dimensions or predictor structure")
    broken <- fit; broken@U <- fit@U[-1,,drop=FALSE]
    expect_error(lfmm2.test(broken,Y,env), "latent factors")
    # Other loci remain allowed; numeric equality / identity of individuals is
    # intentionally a user responsibility, not an unimplemented hash contract.
    expect_length(lfmm2.test(fit,Y[,1:3],env)$pvalues, 6L)
    env$x <- env$x + 0.001
    expect_silent(lfmm2.test(fit,Y,env))
})

test_that("invalid formulas are rejected before any user expression executes", {
    env <- formula_env(); side_effect <- FALSE
    invalid <- list(y ~ x, ~ ., ~ 1, ~ 0, ~ x:y, ~ x*y, ~ x/y, ~ (x+y)^2,
                    ~ scale(x), ~ poly(x,2), ~ cbind(x,y), ~ I(x*y),
                    ~ I(x^0), ~ I(x^-1), ~ I(x^1.5), ~ I(x^y),
                    ~ I(x^(1+1)), ~ log(x,base=2), ~ mean(x), ~ x-x,
                    ~ {side_effect <<- TRUE; x}, ~ get("x"), ~ base::log(w),
                    ~ x + x, ~ sqrt(log(w)))
    for (f in invalid) expect_error(LEA:::.lfmm2_prepare_env(env,f))
    expect_false(side_effect)
    expect_error(LEA:::.lfmm2_prepare_env(env,"~ x"), "Unsupported")
    expect_error(LEA:::.lfmm2_prepare_env(env,~ missing), "Missing formula variables")
    expect_error(LEA:::.lfmm2_prepare_env(as.matrix(env),~ x), "data.frame")
    expect_error(LEA:::.lfmm2_prepare_env("not.env",~ x), "data.frame")
    for (column in list(factor(c("a","b")), c("a","b"), c(TRUE,FALSE),
                        as.Date(c("2000-01-01","2001-01-01")))) {
        bad <- env; bad$bad <- rep(column,length.out=nrow(env))
        expect_error(LEA:::.lfmm2_prepare_env(bad), "only ordinary numeric")
        expect_error(LEA:::.lfmm2_prepare_env(bad,~ x), "only ordinary numeric")
    }
    for (v in c(NA_real_,NaN,Inf,-Inf)) {
        bad <- env; bad$x[1] <- v
        expect_error(LEA:::.lfmm2_prepare_env(bad,~ x), "Non-finite")
    }
    bad <- env; bad$w[1] <- 0
    expect_error(LEA:::.lfmm2_prepare_env(bad,~ log(w)), "finite numeric")
    bad$w[1] <- -1
    expect_error(LEA:::.lfmm2_prepare_env(bad,~ sqrt(w)), "finite numeric")
    bad <- env; bad$x <- 1
    expect_error(LEA:::.lfmm2_prepare_env(bad,~ x), "Constant")
    bad <- env; bad$copy <- 2 * bad$x
    expect_error(LEA:::.lfmm2_prepare_env(bad,~ x + copy), "Linearly dependent")
    bad <- env; bad$x <- rep(0:1,40)
    expect_error(LEA:::.lfmm2_prepare_env(bad,~ x + I(x^2)), "Linearly dependent")
    bad <- env; bad$copy <- bad$x + 1e-7 * bad$z
    expect_silent(LEA:::.lfmm2_prepare_env(bad,~ x + copy))
    names(bad)[2] <- "x"
    expect_error(LEA:::.lfmm2_prepare_env(bad,~ x), "unique")
})

test_that("new S4 objects serialize metadata without changing their five slots", {
    env <- formula_env(); Y <- lfmm_fixtures$binary$Y
    for (f in list(NULL, ~ x + log(w))) {
        fit <- lfmm2(Y,env,2L,effect.sizes=TRUE,formula=f)
        expect_identical(slotNames(fit),c("K","lambda","B","U","V"))
        expect_true(validObject(fit))
        path <- tempfile(); on.exit(unlink(path),add=TRUE)
        saveRDS(fit,path,version=2); copy <- readRDS(path)
        expect_identical(copy,fit)
        save(fit,file=path,version=2); e <- new.env();load(path,e)
        expect_identical(e$fit,fit)
        expect_identical(lfmm2.test(copy,Y,env),lfmm2.test(fit,Y,env))
    }
    bad <- fit; attr(bad,"env.spec") <- list(version=9L)
    expect_error(lfmm2.test(bad,Y,env), "metadata")
})

test_that("real pre-interface RDS and RData S4 fixtures remain valid and usable", {
    root <- test_path("..","lfmm2-s4")
    expected <- readRDS(file.path(root,"expected.rds"))
    e <- new.env();load(file.path(root,"historical.RData"),e)
    objects <- c(e$objects, setNames(lapply(c("no_effect","effect","single"), function(n)
        readRDS(file.path(root,paste0(n,".rds")))), c("no_effect","effect","single")))
    for (index in seq_along(objects)) {
        i <- ((index - 1L) %% 3L) + 1L
        fit <- objects[[index]]
        expect_null(attr(fit,"env.spec",exact=TRUE))
        expect_true(is(fit,"lfmm2Class")); expect_true(validObject(fit))
        for (n in names(expected$slots[[i]])) expect_identical(slot(fit,n),expected$slots[[i]][[n]])
        X <- if(i==3) expected$X[,1,drop=FALSE] else expected$X
        actual <- lapply(c("linear","full","glm"),function(m)
            lfmm2.test(fit,expected$Y,X,full=m=="full",linear=m!="glm"))
        # The fixture was calculated on macOS/R 4.6.1 at 7443970.
        # Preserve attributes exactly, with the existing cross-BLAS numerical
        # tolerances for frozen inference; serialized statistical slots above
        # remain bitwise identical. Same-process matrix/formula checks are exact.
        for (k in seq_along(actual))
            expect_result_close(actual[[k]], expected$results[[i]][[k]], glm=k==3L)
        copy <- unserialize(serialize(fit,NULL,version=2))
        expect_true(validObject(copy)); expect_identical(copy,fit)
    }
})

test_that("formula NULL preserves numeric data-frame and positional fit paths", {
    env <- formula_env(); Y <- lfmm_fixtures$binary$Y
    fit <- lfmm2(Y,env,2L,0.1,TRUE)
    manual <- lfmm2(Y,as.matrix(env),2L,0.1,TRUE)
    for (s in c("K","lambda","B","U","V")) expect_identical(slot(fit,s),slot(manual,s))
    expect_identical(lfmm2.test(fit,Y,env),lfmm2.test(manual,Y,as.matrix(env)))
    expect_null(attr(fit,"env.spec")$formula)
    # The old path deliberately remains available without formula metadata.
    old <- fit; attr(old,"env.spec") <- NULL
    expect_true(validObject(old))
    expect_identical(lfmm2.test(old,Y,env),lfmm2.test(fit,Y,env))
})

test_that("formula language and metadata do not capture executable user bindings", {
    env <- formula_env(); Y <- lfmm_fixtures$binary$Y
    f <- local({ log <- function(...) stop("unexpected user function"); ~ log(w) })
    fit <- lfmm2(Y,env,2L,formula=f)
    X <- matrix(base::log(env$w),ncol=1,dimnames=list(rownames(env),"log(w)"))
    manual <- lfmm2(Y,X,2L)
    expect_identical(fit@U,manual@U)
    expect_identical(lfmm2.test(fit,Y,env),lfmm2.test(manual,Y,X))
    funny <- env; names(funny)[1] <- "with space"
    X <- LEA:::.lfmm2_prepare_env(funny,~ `with space` + y)$X
    expect_identical(colnames(X),c("`with space`","y"))
    fit <- lfmm2(Y,funny,2L,formula=~ `with space` + y)
    expect_identical(lfmm2.test(fit,Y,funny),lfmm2.test(lfmm2(Y,X,2L),Y,X))
    bad <- env; bad$x[1] <- 1e308
    expect_error(LEA:::.lfmm2_prepare_env(bad,~ I(x^3)),"finite numeric")
})

test_that("diagnostic rank is stable under levels and units without modifying raw X", {
    for (x in list(1e16 + seq_len(100), 1 + seq_len(100)*1e-16)) {
        prepared <- LEA:::.lfmm2_prepare_env(data.frame(x=x),~x)$X
        expect_identical(as.vector(prepared), x)
        Y <- lfmm_fixtures$multi$Y[rep(seq_len(80),length.out=100),1:20]
        actual <- lfmm2(Y,data.frame(x=x),2L,formula=~x)
        manual <- lfmm2(Y,prepared,2L)
        for (s in c("U","V","B")) expect_identical(slot(actual,s),slot(manual,s))
        expect_identical(lfmm2.test(actual,Y,data.frame(x=x),full=TRUE),
                         lfmm2.test(manual,Y,prepared,full=TRUE))
    }
    env <- formula_env()
    env$correlated <- env$x + 1e-6*env$z
    expect_silent(LEA:::.lfmm2_prepare_env(env,~x+correlated))
    tiny <- data.frame(a=env$x*1e-100,b=env$y*1e100)
    reference <- as.matrix(tiny); rownames(reference) <- rownames(tiny)
    expect_identical(LEA:::.lfmm2_prepare_env(tiny,~a+b)$X,reference)
    env$exact <- 2*env$x
    expect_error(LEA:::.lfmm2_prepare_env(env,~x+exact),"Linearly dependent")
    env$constant <- 3
    expect_error(LEA:::.lfmm2_prepare_env(env,~constant),"Constant")
    # Subtracting two finite extremes may overflow; the diagnostic handles it.
    extreme <- data.frame(x=rep(c(-1e308,1e308),40))
    reference <- as.matrix(extreme); rownames(reference) <- rownames(extreme)
    expect_identical(LEA:::.lfmm2_prepare_env(extreme,~x)$X,reference)
})

test_that("named transformation arguments fail explicitly before evaluation", {
    env <- formula_env()
    for (f in list(~log(base=x),~sqrt(foo=x),~I(foo=x^2),
                   ~log(x=x),~I(`^`(e1=x,e2=2))))
        expect_error(LEA:::.lfmm2_prepare_env(env,f),"Named arguments")
})

test_that("effective model rank never substitutes a latent statistic for environment", {
    env <- formula_env(); Y <- lfmm_fixtures$binary$Y
    env$close <- env$x + 1e-10*env$z
    fit <- lfmm2(Y,env,2L,formula=~x+close)
    expect_error(lfmm2.test(fit,Y,env,genomic.control=FALSE),
                 "Non-estimable.*close")
    expect_error(lfmm2.test(fit,Y,env,full=TRUE,genomic.control=FALSE),
                 "Non-estimable.*close")
    # Gaussian GLMs use a different fitted QR tolerance and may still retain
    # both terms: rejection follows the actual fit, not the lm default.
    X <- LEA:::.lfmm2_prepare_env(env,~x+close)$X
    raw <- glm(Y[,1]~.,data=data.frame(X,fit@U),family=gaussian())
    if (raw$rank==length(raw$coefficients))
        expect_identical(lfmm2.test(fit,Y,env,linear=FALSE,family=gaussian(),genomic.control=FALSE),
                         lfmm2.test(lfmm2(Y,X,2L),Y,X,linear=FALSE,family=gaussian(),genomic.control=FALSE))
    # A tighter dependence falls below the GLM tolerance too.
    env$close <- env$x + 1e-13*env$z
    fit <- lfmm2(Y,env,2L,formula=~x+close)
    expect_error(lfmm2.test(fit,Y,env,linear=FALSE,family=gaussian(),genomic.control=FALSE),
                 "Non-estimable.*close")
})

test_that("confounding with U is rejected but redundant nuisance factors alone are allowed", {
    env <- formula_env(); Y <- lfmm_fixtures$binary$Y
    fit <- lfmm2(Y,env,2L,formula=~x+y)
    fit@U[,1] <- env$x
    for (mode in c("linear","full","glm"))
        expect_error(lfmm2.test(fit,Y,env,full=mode=="full",linear=mode!="glm",
                               family=gaussian(),genomic.control=FALSE),"Non-estimable.*x")
    # Legacy behavior is explicitly preserved for otherwise identical objects.
    legacy <- fit; attr(legacy,"env.spec") <- NULL
    X <- LEA:::.lfmm2_prepare_env(env,~x+y)$X
    for (mode in c("linear","full","glm"))
        expect_result_close(lfmm2.test(legacy,Y,X,full=mode=="full",linear=mode!="glm",
                                       family=gaussian(),genomic.control=FALSE),
                            historical_test(legacy,Y,X,full=mode=="full",linear=mode!="glm",
                                            family=gaussian(),genomic.control=FALSE),glm=mode=="glm")
    fit <- lfmm2(Y,env,2L,formula=~x+y)
    fit@U[,2] <- fit@U[,1]
    manual <- fit; attr(manual,"env.spec") <- NULL
    for (mode in c("linear","full","glm"))
        expect_identical(lfmm2.test(fit,Y,env,full=mode=="full",linear=mode!="glm",
                                    family=gaussian(),genomic.control=FALSE),
                         lfmm2.test(manual,Y,X,full=mode=="full",linear=mode!="glm",
                                    family=gaussian(),genomic.control=FALSE))
})

test_that("GLM estimability is checked again when a later locus has different weights", {
    env <- formula_env(); fit <- lfmm2(lfmm_fixtures$binary$Y,env,2L,formula=~x+y)
    Y <- cbind(seq(-2,2,length.out=80),seq(1,3,length.out=80))
    family <- gaussian()
    # An initializer may select response-dependent prior weights. The first
    # locus retains all rows; the second has only two informative observations.
    family$initialize <- expression({
        n <- rep.int(1,nobs)
        if (all(y>0)) weights <- as.numeric(seq_len(nobs)<=2L)
        mustart <- y
    })
    expect_false(LEA:::.lfmm2_glm_reusable_family(family))
    expect_error(lfmm2.test(fit,Y,env,linear=FALSE,family=family,genomic.control=FALSE),
                 "Non-estimable.*y")
})

test_that("full-rank correlated formula designs retain reference inference", {
    env <- formula_env(); Y <- lfmm_fixtures$binary$Y
    env$close <- env$x + 1e-6*env$z
    fit <- lfmm2(Y,env,2L,formula=~x+close)
    X <- LEA:::.lfmm2_prepare_env(env,~x+close)$X
    manual <- lfmm2(Y,X,2L)
    for (mode in c("linear","full","glm")) for (gc in c(FALSE,TRUE))
        expect_identical(lfmm2.test(fit,Y,env,full=mode=="full",linear=mode!="glm",
                                    genomic.control=gc),
                         lfmm2.test(manual,Y,X,full=mode=="full",linear=mode!="glm",
                                    genomic.control=gc))
})

test_that("global estimability follows the retained nuisance and residual models", {
    env <- formula_env(); Y <- lfmm_fixtures$binary$Y
    fit <- lfmm2(Y,env,2L,formula=~x+y)
    fit@U[,1] <- env$x + 1e-8*env$z
    manual <- fit; attr(manual,"env.spec") <- NULL
    X <- LEA:::.lfmm2_prepare_env(env,~x+y)$X
    expect_identical(lfmm2.test(fit,Y,env,full=TRUE),
                     lfmm2.test(manual,Y,X,full=TRUE))
    # The first stage drops the almost duplicated second nuisance column.
    # Do not restore it in the diagnostic and spuriously absorb x.
    fit <- lfmm2(Y,env,2L,formula=~x+y)
    fit@U[,2] <- fit@U[,1] + 1e-10*env$x
    manual <- fit; attr(manual,"env.spec") <- NULL
    expect_identical(lfmm2.test(fit,Y,env,full=TRUE),
                     lfmm2.test(manual,Y,X,full=TRUE))
})

test_that("env.spec rejects malformed keys, versions, field types and dimensions", {
    env <- formula_env(); Y <- lfmm_fixtures$binary$Y
    fit <- lfmm2(Y,env,2L,formula=~x+I(y^2))
    valid <- attr(fit,"env.spec",exact=TRUE)
    replace <- function(name,value) {s<-valid;s[name]<-list(value);s}
    invalid <- list(NULL_formula=replace("formula",NULL),
        unknown_version=replace("version",2L),double_version=replace("version",1),
        missing_version=valid[names(valid)!="version"],
        partial_version=setNames(valid,sub("^version$","version_future",names(valid))),
        duplicate_key=c(valid,list(version=1L)),unexpected_key=c(valid,list(extra=1L)),
        missing_names=unname(valid),duplicate_other=c(valid,list(mode="formula")),
        bad_name=setNames(valid,c(NA_character_,names(valid)[-1])),
        mode_type=replace("mode",1L),mode_NA=replace("mode",NA_character_),
        mode_multiple=replace("mode",c("formula","matrix")),mode_unknown=replace("mode","unknown"),
        n_double=replace("n",as.double(valid$n)),n_NA=replace("n",NA_integer_),
        n_zero=replace("n",0L),n_multiple=replace("n",c(80L,80L)),
        n_named=replace("n",c(n=80L)),n_list=replace("n",list(80L)),
        d_negative=replace("d",-1L),d_zero=replace("d",0L),
        d_wide=replace("d",valid$n+1L),d_wrong=replace("d",1L),
        predictors_type=replace("predictors",1:2),predictors_NULL=replace("predictors",NULL),
        predictors_length=replace("predictors","x"),predictors_NA=replace("predictors",c("x",NA)),
        predictors_duplicates=replace("predictors",c("x","x")),
        predictors_empty=replace("predictors",c("x","")),
        predictors_order=replace("predictors",rev(valid$predictors)),
        predictors_matrix=replace("predictors",matrix(valid$predictors,nrow=1)),
        formula_type=replace("formula",~x+I(y^2)),formula_empty=replace("formula",""),
        formula_NA=replace("formula",NA_character_),
        formula_multiple=replace("formula",c("~x","~y")),
        formula_syntax=replace("formula","~ x +"),formula_response=replace("formula","y ~ x"),
        formula_unsupported=replace("formula","~ scale(x)"),
        formula_order=replace("formula","~ I(y^2) + x"),
        formula_array=replace("formula",matrix("~ x")))
    for (name in names(invalid)) {
        bad <- fit; attr(bad,"env.spec") <- invalid[[name]]
        expect_error(LEA:::.lfmm2_env_spec(bad),"metadata",info=name)
        expect_error(lfmm2.test(bad,Y,env),"metadata",info=name)
    }
    # These are present attributes, not the historical absent-attribute case.
    for (value in list(1L,list(),structure(valid,class="custom"),as.pairlist(valid))) {
        bad <- fit;attr(bad,"env.spec") <- value
        expect_error(LEA:::.lfmm2_env_spec(bad),"metadata")
    }
    s <- valid; s$n <- s$d
    bad <- fit;attr(bad,"env.spec") <- s
    expect_error(LEA:::.lfmm2_env_spec(bad),"metadata")
    # Parsing a formula field must never evaluate a non-formula expression.
    old <- getOption("LEA_metadata_eval");on.exit(options(LEA_metadata_eval=old),add=TRUE)
    options(LEA_metadata_eval=NULL)
    s <- replace("formula","{options(LEA_metadata_eval=TRUE); ~x}")
    bad <- fit;attr(bad,"env.spec") <- s
    expect_error(LEA:::.lfmm2_env_spec(bad),"metadata")
    expect_null(getOption("LEA_metadata_eval"))
})

test_that("env.spec preserves historical objects and valid input modes and matrix names", {
    env <- formula_env();Y <- lfmm_fixtures$binary$Y
    fit <- lfmm2(Y,env,2L,formula=~x+I(y^2))
    spec <- attr(fit,"env.spec",exact=TRUE)
    expect_identical(LEA:::.lfmm2_env_spec(fit),spec)
    reordered <- fit;attr(reordered,"env.spec") <- rev(spec)
    expect_identical(LEA:::.lfmm2_env_spec(reordered),rev(spec))
    expect_identical(LEA:::.lfmm2_env_spec(unserialize(serialize(fit,NULL,version=2))),spec)
    old <- readRDS(test_path("..","lfmm2-s4","no_effect.rds"))
    expect_null(LEA:::.lfmm2_env_spec(old));expect_true(validObject(old))
    for (labels in list(NULL,c("duplicate","duplicate"),c("",NA_character_))) {
        X <- as.matrix(env[,c("x","y")]);colnames(X) <- labels
        object <- lfmm2(Y,X,2L)
        expect_identical(LEA:::.lfmm2_env_spec(object)$predictors,labels)
        expect_null(LEA:::.lfmm2_env_spec(object)$formula)
        if (is.null(labels) || all(!is.na(labels) & nzchar(labels))) {
            historical <- object;attr(historical,"env.spec") <- NULL
            expect_identical(lfmm2.test(object,Y,X),lfmm2.test(historical,Y,X))
        }
    }
    for (E in list(as.matrix(env),as.vector(env$x),env)) {
        object <- lfmm2(Y,E,2L)
        expect_identical(LEA:::.lfmm2_env_spec(object),attr(object,"env.spec"))
    }
    file <- tempfile(fileext=".env");on.exit(unlink(file),add=TRUE)
    write.table(as.matrix(env),file,row.names=FALSE,col.names=FALSE)
    object <- lfmm2(Y,file,2L)
    expect_identical(LEA:::.lfmm2_env_spec(object)$mode,"file")
    # No new equality restriction for no-formula metadata: another valid env
    # remains accepted exactly as on the historical association path.
    historical <- object;attr(historical,"env.spec") <- NULL
    expect_identical(lfmm2.test(object,Y,env$x),lfmm2.test(historical,Y,env$x))
    for (mode in c("matrix","data.frame","file")) {
        bad <- fit;s<-spec;s$mode<-mode;attr(bad,"env.spec")<-s
        expect_error(LEA:::.lfmm2_env_spec(bad),"metadata")
    }
    # d=n is valid metadata without a formula, not on a centered formula path.
    X <- matrix(c(0,1,1,0),2);object <- lfmm2(matrix(1:6,2),X,1L)
    expect_identical(LEA:::.lfmm2_env_spec(object)$d,2L)
})
