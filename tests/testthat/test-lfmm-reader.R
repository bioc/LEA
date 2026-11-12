reader_capture <- function(fun, input) {
    warnings <- character()
    value <- tryCatch(withCallingHandlers(fun(input), warning = function(w) {
        warnings <<- c(warnings, gsub(input, "<INPUT>", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    }), error = identity)
    list(ok = !inherits(value, "error"), value = value, warnings = warnings)
}

reader_fixture <- function(name) {
    input <- tempfile(fileext = ".lfmm")
    bytes <- readRDS(test_path("..", "lfmm-reader", "fixtures.rds"))[[paste0(name, ".txt")]]
    stopifnot(is.raw(bytes))
    writeBin(bytes, input)
    input
}

test_that("published reader fixtures preserve exact values, types, names and warnings", {
    expected <- readRDS(test_path("..", "lfmm-reader", "expected.rds"))
    for (name in names(expected)) {
        input <- reader_fixture(name)
        actual <- reader_capture(read.lfmm, input)
        expect_identical(actual$ok, expected[[name]]$ok, info = name)
        expect_identical(actual$warnings, expected[[name]]$warnings, info = name)
        if (actual$ok) expect_identical(actual$value, expected[[name]]$value, info = name)
        else {
            legacy <- reader_capture(function(p) as.matrix(read.table(p)), input)
            expect_identical(conditionMessage(actual$value), conditionMessage(legacy$value), info = name)
        }
        unlink(input)
    }
})

test_that("regular numeric inputs use the checked fast path", {
    for (name in c("regular", "tabs", "mixed", "crlf", "blanks", "no-final-newline",
                   "single-column", "dosage", "integral-dosage", "integer-spelling",
                   "missing", "all-NA", "NA-column", "limits", "beyond-integer",
                   "hex", "exponents")) {
        input <- reader_fixture(name)
        actual <- LEA:::.lfmm_read_numeric(input)
        expect_true(is.matrix(actual), info = name)
        expect_identical(actual, suppressWarnings(as.matrix(read.table(input))), info = name)
        unlink(input)
    }
})

test_that("ambiguous legacy syntax and malformed records never get blindly reshaped", {
    for (name in c("single-row", "short-no-newline", "comment", "quoted", "logical",
                   "logical-numeric", "bad-token", "complex", "header-rowids", "nonfinite",
                   "ragged", "balanced-ragged", "extra-column", "empty", "comments-only")) {
        input <- reader_fixture(name)
        expect_null(LEA:::.lfmm_read_numeric(input), info = name)
        unlink(input)
    }
    for (name in c("ragged", "balanced-ragged", "extra-column", "empty", "comments-only")) {
        input <- reader_fixture(name)
        expect_error(read.lfmm(input))
        unlink(input)
    }
})

test_that("block boundaries, late dosage and integer limits preserve complete records", {
    for (late in c("2", "1.0", "2147483648", "NA")) {
        input <- tempfile(fileext = ".lfmm")
        lines <- rep("0 1 2", 600)
        lines[599] <- paste("0 1", late)
        writeLines(c("", lines[1:256], " \t", lines[257:600]), input)
        expect_identical(read.lfmm(input), as.matrix(read.table(input)))
        expect_identical(LEA:::.lfmm_read_numeric(input), as.matrix(read.table(input)))
        unlink(input)
    }
})

test_that("numeric token spellings match the original conversion rules", {
    set.seed(8529)
    tokens <- c("0", "1", "2", "-9", "NA", "0.0", "1e0", "0x10", "1e-320",
                "+001", "2147483648", "-2147483648", "Inf", "NaN")
    for (i in 1:20) {
        input <- tempfile(fileext = ".lfmm")
        x <- matrix(sample(tokens, 48, replace = TRUE), 12, 4)
        writeLines(apply(x, 1, paste, collapse = " "), input)
        expect_identical(read.lfmm(input), as.matrix(read.table(input)))
        unlink(input)
    }
})

test_that("extension validation and repeated fallback errors retain their contracts", {
    expect_error(read.lfmm(), "missing")
    expect_error(read.lfmm(1), "character")
    expect_error(read.lfmm("input.env"), "extension")
    before <- nrow(showConnections(all = TRUE))
    input <- reader_fixture("balanced-ragged")
    for (i in 1:30) expect_error(read.lfmm(input))
    unlink(input)
    expect_identical(nrow(showConnections(all = TRUE)), before)
})

test_that("non-ASCII and control bytes are delegated without extra diagnostics", {
    for (byte in c(1L, 11L, 12L, 127L, 128L, 255L)) {
        input <- tempfile(fileext = ".lfmm")
        line <- c(charToRaw("0 "), as.raw(byte), charToRaw(" 1\n"))
        writeBin(rep(line, 8), input)
        actual <- reader_capture(read.lfmm, input)
        legacy <- reader_capture(function(p) as.matrix(read.table(p)), input)
        expect_identical(actual$ok, legacy$ok)
        expect_identical(actual$warnings, legacy$warnings)
        if (actual$ok) expect_identical(actual$value, legacy$value)
        else expect_identical(conditionMessage(actual$value), conditionMessage(legacy$value))
        unlink(input)
    }
})

test_that("permissive scan spellings cannot silently change legacy token types", {
    tokens <- c("+NA", "-NA", "NAN", "NANAN", "NANANA", "A", "NA1", "+NaN", "-NaN",
                "infinity", "Infinity", "-Infinity", "nan", "1D3", "1d3",
                paste(rep("9", 400), collapse = ""))
    for (token in tokens) {
        input <- tempfile(fileext = ".lfmm")
        writeLines(rep(paste(token, "0 1"), 8), input)
        actual <- reader_capture(read.lfmm, input)
        legacy <- reader_capture(function(p) as.matrix(read.table(p)), input)
        expect_identical(actual$ok, legacy$ok, info = token)
        expect_identical(actual$warnings, legacy$warnings, info = token)
        if (actual$ok) expect_identical(actual$value, legacy$value, info = token)
        else expect_identical(conditionMessage(actual$value), conditionMessage(legacy$value), info = token)
        unlink(input)
    }
})
