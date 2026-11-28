#' @noRd
#' @export
rgraphtool_max_threads <- function() {
  as.integer(cpp_max_threads())
}

#' @noRd
#' @export
rgraphtool_num_threads <- function() {
  opt <- getOption("rgraphtool.nthreads", NULL)
  if (!is.null(opt)) {
    if (is.character(opt) && length(opt) == 1L && identical(tolower(opt), "auto")) {
      return(rgraphtool_max_threads())
    }
    if (is.numeric(opt) && length(opt) == 1L && !is.na(opt) &&
        is.finite(opt) && opt >= 1 && opt == floor(opt)) {
      return(as.integer(min(opt, 256)))
    }
    stop("option rgraphtool.nthreads must be a positive integer or \"auto\"",
         call. = FALSE)
  }
  as.integer(cpp_get_num_threads())
}

#' @noRd
#' @export
rgraphtool_set_num_threads <- function(nthreads) {
  if (is.character(nthreads) && length(nthreads) == 1L &&
      identical(tolower(nthreads), "auto")) {
    n <- rgraphtool_max_threads()
  } else if (is.numeric(nthreads) && length(nthreads) == 1L && !is.na(nthreads) &&
             is.finite(nthreads) && nthreads >= 1 && nthreads == floor(nthreads)) {
    n <- as.integer(min(nthreads, 256))
  } else {
    stop("nthreads must be a positive integer or \"auto\"", call. = FALSE)
  }
  invisible(as.integer(cpp_set_num_threads(n)))
}

resolve_nthreads <- function(nthreads = NULL) {
  if (is.null(nthreads)) {
    return(rgraphtool_num_threads())
  }
  if (is.character(nthreads) && length(nthreads) == 1L &&
      identical(tolower(nthreads), "auto")) {
    return(rgraphtool_max_threads())
  }
  if (!is.numeric(nthreads) || length(nthreads) != 1L || is.na(nthreads) ||
      !is.finite(nthreads) || nthreads < 1 || nthreads != floor(nthreads)) {
    stop("nthreads must be a positive integer, \"auto\", or NULL", call. = FALSE)
  }
  as.integer(min(nthreads, 256))
}

# Number of worker processes for independent starts: NULL uses
# getOption("rgraphtool.cores", 1L). Forked workers are unavailable on Windows.
resolve_cores <- function(cores = NULL) {
  if (is.null(cores)) cores <- getOption("rgraphtool.cores", 1L)
  if (!is.numeric(cores) || length(cores) != 1L || is.na(cores) || cores < 1 ||
      cores != floor(cores)) {
    stop("cores must be NULL or a positive integer", call. = FALSE)
  }
  cores <- as.integer(cores)
  if (cores > 1L && .Platform$OS.type == "windows") {
    warning("cores > 1 needs forked processes, which are unavailable on Windows; ",
            "running starts sequentially", call. = FALSE)
    cores <- 1L
  }
  cores
}

# Replace every (serialized, pointer-less) rgraphtool_graph inside `x` with
# the live `graph`.
reattach_graph <- function(x, graph) {
  if (inherits(x, "rgraphtool_graph")) return(graph)
  if (is.list(x) && !is.data.frame(x)) {
    at <- attributes(x)
    x <- lapply(x, reattach_graph, graph = graph)
    attributes(x) <- at
  }
  x
}

#' @noRd
#' @export
parallel_starts <- function(n, fun, cores = NULL, graph = NULL) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1 || n != floor(n)) {
    stop("n must be a positive integer", call. = FALSE)
  }
  n <- as.integer(n)
  if (!is.function(fun)) stop("fun must be a function of the run index", call. = FALSE)
  cores <- resolve_cores(cores)
  if (n == 1L && cores == 1L) return(list(fun(1L)))
  # Per-run seeds drawn up front: results are identical for any `cores`.
  seeds <- sample.int(.Machine$integer.max, n)
  after_seeds <- get(".Random.seed", envir = globalenv())
  run_one <- function(i) {
    set.seed(seeds[i])
    fun(i)
  }
  if (cores == 1L) {
    res <- lapply(seq_len(n), run_one)
    assign(".Random.seed", after_seeds, envir = globalenv())
    return(res)
  }
  child <- function(i) {
    # The parent's persistent thread pool does not exist in a forked child:
    # workers must run single-threaded.
    options(rgraphtool.nthreads = 1L)
    cpp_set_num_threads(1L)
    run_one(i)
  }
  res <- withCallingHandlers(
    parallel::mclapply(seq_len(n), child, mc.cores = min(cores, n),
                       mc.preschedule = FALSE, mc.set.seed = FALSE),
    # failures are reported below as an error; drop mclapply's own warning
    warning = function(w) {
      if (grepl("function calls? resulted in an error", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  failed <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(failed)) {
    msg <- res[[which(failed)[1L]]]
    stop("a parallel start failed: ",
         if (inherits(msg, "try-error")) conditionMessage(attr(msg, "condition")) else "no result",
         call. = FALSE)
  }
  if (!is.null(graph)) res <- lapply(res, reattach_graph, graph = graph)
  res
}
