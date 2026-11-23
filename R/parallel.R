rgraphtool_max_threads <- function() {
  as.integer(cpp_max_threads())
}

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
