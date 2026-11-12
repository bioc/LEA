read.lfmm <- function(input.file) {

        # test arguments
        if(missing(input.file))
                stop("'input.file' argument is missing.")
        else if (!is.character(input.file))
                stop("'input.file' argument has to be of type character.")
    # check extension 
    test_extension(input.file, "lfmm")

    numeric <- .lfmm_read_numeric(input.file)
    if (!is.null(numeric)) return(numeric)
    return(as.matrix(read.table(input.file)))
}

# Return NULL for legacy syntax. Structural validation precedes numeric scanning;
# the legacy parser remains the authority for ambiguous or malformed input.
.lfmm_read_numeric <- function(input.file) {
    if (!utils::file_test("-f", input.file)) return(NULL)
    signature <- function() file.info(input.file)[, c("size", "mtime", "ctime"), drop = FALSE]
    before <- signature()
    warned <- FALSE
    widths <- tryCatch(withCallingHandlers(
        utils::count.fields(input.file, quote = "\"'", comment.char = "#",
                     blank.lines.skip = TRUE),
        warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }),
        error = function(e) NULL)
    if (warned || !length(widths) || anyNA(widths) || any(widths < 1L) ||
        any(widths != widths[1L])) return(NULL)
    n <- length(widths)
    p <- widths[1L]
    # Keep read.table's short-file diagnostics, including incomplete final lines.
    if (n <= 5L) return(NULL)
    if (n > .Machine$integer.max || p > .Machine$integer.max ||
        !is.finite(as.double(n) * p)) return(NULL)

    con <- file(input.file, "rt")
    on.exit(close(con))
    output <- matrix(NA_real_, nrow = n, ncol = p)
    chunk <- max(1L, min(256L, floor(1e6 / p)))
    done <- 0L
    integer_tokens <- TRUE
    digits <- FALSE
    repeat {
        lines <- tryCatch(withCallingHandlers(readLines(con, n = chunk, warn = FALSE),
            warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }),
            error = function(e) NULL)
        if (is.null(lines) || warned) return(NULL)
        if (!length(lines)) break
        # Non-ASCII/control bytes belong to the legacy parser, without first
        # applying locale-sensitive regular expressions to malformed encodings.
        if (any(grepl("[^\t -~]", lines, useBytes = TRUE))) return(NULL)
        if (any(grepl("[#\"']", lines))) return(NULL)
        # scan() accepts nonfinite spellings (for example NAN) which legacy
        # type.convert() can retain as text. Only canonical NA stays on this path.
        if (any(grepl("[NnIi]", lines)) &&
            any(grepl("[NnIi]", gsub("NA", "", lines, fixed = TRUE)))) return(NULL)
        lines <- lines[!grepl("^[[:space:]]*$", lines)]
        rows <- length(lines)
        if (!rows) next
        if (done > n - rows) return(NULL)
        # Lexical integer detection preserves 1.0/1e0 as double, even when their
        # numeric values are integral. scan() validates the actual token grammar.
        integer_tokens <- integer_tokens && !any(grepl("[^0-9+NA \t-]", lines))
        block_digits <- any(grepl("[0-9]", lines))
        digits <- digits || block_digits
        text <- textConnection(lines)
        values <- tryCatch(withCallingHandlers(
            scan(text, what = double(), nmax = as.double(rows) * p + 1,
                 quiet = TRUE),
            warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }),
            error = function(e) NULL, finally = close(text))
        if (is.null(values) || warned || length(values) != as.double(rows) * p)
            return(NULL)
        if (integer_tokens && block_digits &&
            (min(values, na.rm = TRUE) < -.Machine$integer.max ||
             max(values, na.rm = TRUE) > .Machine$integer.max)) integer_tokens <- FALSE
        output[seq.int(done + 1L, done + rows), ] <-
            matrix(values, nrow = rows, ncol = p, byrow = TRUE)
        done <- done + rows
    }
    if (done != n || !identical(before, signature())) return(NULL)
    if (integer_tokens) storage.mode(output) <- if (digits) "integer" else "logical"
    colnames(output) <- paste0("V", seq_len(p))
    output
}
