# R interface for C++ extension algorithms, header discovery, and user
# Rcpp::sourceCpp() compilation helpers.

#' @noRd
#' @export
kcore_decomposition <- function(graph, deg = c("total", "out", "in"), checked = FALSE) {
  graph_pointer(graph)
  deg <- match.arg(deg)
  if (!is.logical(checked) || length(checked) != 1L || is.na(checked)) {
    stop("checked must be TRUE or FALSE", call. = FALSE)
  }
  if (!isTRUE(graph$directed) && deg != "total") {
    stop("Undirected graphs support only deg = 'total'", call. = FALSE)
  }
  cpp_kcore_decomposition(graph, deg_mode = deg, checked = checked)
}

#' @noRd
#' @export
vertex_strength <- function(graph, property, mode = c("total", "out", "in"), checked = FALSE) {
  graph_pointer(graph)
  mode <- match.arg(mode)
  if (!is.character(property) || length(property) != 1L || !nzchar(property)) {
    stop("property must be a non-empty edge property name", call. = FALSE)
  }
  if (!is.logical(checked) || length(checked) != 1L || is.na(checked)) {
    stop("checked must be TRUE or FALSE", call. = FALSE)
  }
  if (!isTRUE(graph$directed) && mode != "total") {
    stop("Undirected graphs support only mode = 'total'", call. = FALSE)
  }
  cpp_vertex_strength(graph, prop_name = property, mode = mode, checked = checked)
}

#' @noRd
#' @export
rgraphtool_include_dir <- function() {
  inc <- system.file("include", package = "rgraphtool")
  if (nzchar(inc) && file.exists(file.path(inc, "rgraphtool.h"))) {
    return(normalizePath(inc, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate rgraphtool's include directory containing rgraphtool.h",
       call. = FALSE)
}

#' @noRd
#' @export
rgraphtool_extension_example <- function(name = c("kcore", "weighted_degree")) {
  name <- match.arg(name)
  fname <- switch(name,
    kcore = "kcore_extension.cpp",
    weighted_degree = "weighted_degree_extension.cpp"
  )
  path <- system.file("examples", "extensions", fname, package = "rgraphtool")
  if (nzchar(path) && file.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  stop(sprintf("Could not locate bundled C++ extension example '%s'", fname),
       call. = FALSE)
}

#' @noRd
#' @export
compile_cpp_extension <- function(file = NULL, code = NULL,
                                  env = parent.frame(),
                                  rebuild = FALSE,
                                  verbose = FALSE, ...) {
  if (is.null(file) == is.null(code)) {
    stop("Supply exactly one of 'file' or 'code'", call. = FALSE)
  }
  if (!is.null(file)) {
    if (!is.character(file) || length(file) != 1L || !file.exists(file)) {
      stop("file must point to an existing C++ source file", call. = FALSE)
    }
    src_lines <- readLines(file, warn = FALSE)
  } else {
    if (!is.character(code) || length(code) < 1L) {
      stop("code must be a character vector of C++ source code", call. = FALSE)
    }
    src_lines <- unlist(strsplit(paste(code, collapse = "\n"), "\n", fixed = TRUE))
  }

  # Strip 'rgraphtool' from Rcpp::depends so sourceCpp works both when
  # rgraphtool is installed and when loaded via pkgload::load_all(), while
  # explicitly passing rgraphtool_include_dir() via PKG_CPPFLAGS.
  src_lines <- vapply(src_lines, function(line) {
    if (grepl("^\\s*//\\s*\\[\\[\\s*Rcpp::depends\\(", line)) {
      line <- gsub("\\brgraphtool\\b\\s*,\\s*", "", line)
      line <- gsub(",\\s*\\brgraphtool\\b", "", line)
      line <- gsub("Rcpp::depends\\(\\s*rgraphtool\\s*\\)", "Rcpp::depends(BH)", line)
    }
    line
  }, character(1L), USE.NAMES = FALSE)

  # Ensure BH dependency is declared if not already present
  has_bh <- any(grepl("Rcpp::depends\\(.*\\bBH\\b", src_lines))
  has_cpp_std <- any(grepl("Rcpp::plugins\\(.*cpp(17|20|23)", src_lines))
  header_lines <- character(0)
  if (!has_bh) header_lines <- c(header_lines, "// [[Rcpp::depends(BH)]]")
  if (!has_cpp_std) header_lines <- c(header_lines, "// [[Rcpp::plugins(cpp23)]]")
  full_code <- paste(c(header_lines, src_lines), collapse = "\n")

  inc_dir <- rgraphtool_include_dir()
  old_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = NA_character_)
  inc_flag <- sprintf("-I\"%s\"", inc_dir)
  new_cppflags <- if (is.na(old_cppflags) || !nzchar(old_cppflags)) {
    inc_flag
  } else {
    paste(old_cppflags, inc_flag)
  }
  Sys.setenv(PKG_CPPFLAGS = new_cppflags)
  on.exit({
    if (is.na(old_cppflags)) {
      Sys.unsetenv("PKG_CPPFLAGS")
    } else {
      Sys.setenv(PKG_CPPFLAGS = old_cppflags)
    }
  }, add = TRUE)

  Rcpp::sourceCpp(code = full_code, env = env, rebuild = rebuild,
                  verbose = verbose, ...)
}
