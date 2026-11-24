read_all_raw <- function(con, chunk_size = 1048576L) {
  chunks <- list()
  i <- 1L
  repeat {
    buf <- readBin(con, what = "raw", n = chunk_size)
    if (!length(buf)) break
    chunks[[i]] <- buf
    i <- i + 1L
  }
  if (!length(chunks)) return(raw(0))
  do.call(c, chunks)
}

detect_gt_compression <- function(file) {
  con <- file(file, open = "rb")
  on.exit(close(con), add = TRUE)
  hdr <- readBin(con, what = "raw", n = 6L)
  if (length(hdr) >= 6L &&
      identical(hdr, as.raw(c(0xe2, 0x9b, 0xbe, 0x20, 0x67, 0x74)))) {
    return("none")
  }
  if (length(hdr) >= 2L && identical(hdr[1:2], as.raw(c(0x1f, 0x8b)))) {
    return("gzip")
  }
  if (length(hdr) >= 3L && identical(hdr[1:3], as.raw(c(0x42, 0x5a, 0x68)))) {
    return("bzip2")
  }
  if (length(hdr) >= 6L &&
      identical(hdr, as.raw(c(0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00)))) {
    return("xz")
  }
  if (length(hdr) >= 4L &&
      identical(hdr[1:4], as.raw(c(0x28, 0xb5, 0x2f, 0xfd)))) {
    return("zstd")
  }
  "none"
}

infer_gt_compression_from_path <- function(file) {
  lf <- tolower(file)
  if (endsWith(lf, ".gz")) return("gzip")
  if (endsWith(lf, ".bz2")) return("bzip2")
  if (endsWith(lf, ".xz")) return("xz")
  if (endsWith(lf, ".zst") || endsWith(lf, ".zstd")) return("zstd")
  "none"
}

#' @noRd
#' @export
read_gt <- function(file) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("file must be a non-empty character path", call. = FALSE)
  }
  path <- path.expand(file)
  if (!file.exists(path) || dir.exists(path)) {
    stop(sprintf("File does not exist: %s", file), call. = FALSE)
  }

  comp <- detect_gt_compression(path)
  parsed <- if (comp == "none") {
    cpp_read_gt_file(path)
  } else if (comp %in% c("gzip", "bzip2", "xz")) {
    con <- switch(comp,
                  gzip = gzfile(path, open = "rb"),
                  bzip2 = bzfile(path, open = "rb"),
                  xz = xzfile(path, open = "rb"))
    on.exit(close(con), add = TRUE)
    cpp_read_gt_raw(read_all_raw(con))
  } else if (comp == "zstd") {
    zstd_bin <- Sys.which("zstd")
    if (!nzchar(zstd_bin)) {
      stop("Reading zstd-compressed .gt files requires the 'zstd' command-line utility in PATH",
           call. = FALSE)
    }
    tmp <- tempfile(fileext = ".gt")
    on.exit(unlink(tmp), add = TRUE)
    status <- system2(zstd_bin, c("-d", "-f", "-q", path, "-o", tmp))
    if (!identical(status, 0L) || !file.exists(tmp)) {
      stop(sprintf("Failed to decompress zstd file: %s", file), call. = FALSE)
    }
    cpp_read_gt_file(tmp)
  }

  attrs <- new_graph_attrs()
  attrs$graph_props <- parsed$graph_props
  attrs$vertex_props <- parsed$vertex_props
  attrs$edge_props <- parsed$edge_props_extra
  attrs$comment <- parsed$comment

  structure(
    list(pointer = parsed$pointer,
         directed = isTRUE(parsed$directed),
         attrs = attrs),
    class = "rgraphtool_graph"
  )
}

#' @noRd
#' @export
write_gt <- function(graph, file, comment = NULL,
                     compression = c("auto", "none", "gzip", "bzip2", "xz", "zstd")) {
  pointer <- graph_pointer(graph)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("file must be a non-empty character path", call. = FALSE)
  }
  compression <- match.arg(compression)
  path <- path.expand(file)
  if (compression == "auto") {
    compression <- infer_gt_compression_from_path(path)
  }

  nv <- graph_num_vertices(graph)
  ne <- graph_num_edges(graph)
  gprops <- if (is.environment(graph$attrs)) graph$attrs$graph_props else list()
  vprops <- if (is.environment(graph$attrs)) graph$attrs$vertex_props else list()
  eprops_extra <- if (is.environment(graph$attrs)) graph$attrs$edge_props else list()

  if (is.null(comment)) {
    n_ep <- length(cpp_edge_property_names(pointer)) + length(eprops_extra)
    comment <- sprintf(
      "rgraphtool binary file version 1 (https://github.com/diogro/rgraphtool), stats: %d vertices, %d edges, %s, %d graph props, %d vertex props, %d edge props",
      as.integer(nv), as.integer(ne),
      if (isTRUE(graph$directed)) "directed" else "undirected",
      length(gprops), length(vprops), n_ep
    )
  } else if (!is.character(comment) || length(comment) != 1L || is.na(comment)) {
    stop("comment must be NULL or a single character string", call. = FALSE)
  }

  if (compression == "none") {
    cpp_write_gt_file(pointer, path, comment, gprops, vprops, eprops_extra)
  } else if (compression %in% c("gzip", "bzip2", "xz")) {
    raw_bytes <- cpp_write_gt_raw(pointer, comment, gprops, vprops, eprops_extra)
    con <- switch(compression,
                  gzip = gzfile(path, open = "wb"),
                  bzip2 = bzfile(path, open = "wb"),
                  xz = xzfile(path, open = "wb"))
    on.exit(close(con), add = TRUE)
    writeBin(raw_bytes, con)
  } else if (compression == "zstd") {
    zstd_bin <- Sys.which("zstd")
    if (!nzchar(zstd_bin)) {
      stop("Writing zstd-compressed .gt files requires the 'zstd' command-line utility in PATH",
           call. = FALSE)
    }
    tmp <- tempfile(fileext = ".gt")
    on.exit(unlink(tmp), add = TRUE)
    cpp_write_gt_file(pointer, tmp, comment, gprops, vprops, eprops_extra)
    status <- system2(zstd_bin, c("-f", "-q", tmp, "-o", path))
    if (!identical(status, 0L)) {
      stop(sprintf("Failed to compress .gt file with zstd: %s", file), call. = FALSE)
    }
  }

  invisible(file)
}

require_igraph <- function(fn_name) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(sprintf("Package 'igraph' is required for %s(); install it with install.packages('igraph')",
                 fn_name), call. = FALSE)
  }
}

#' @noRd
#' @export
as_igraph <- function(graph, include_properties = TRUE) {
  require_igraph("as_igraph")
  graph_pointer(graph)
  if (!is.logical(include_properties) || length(include_properties) != 1L ||
      is.na(include_properties)) {
    stop("include_properties must be TRUE or FALSE", call. = FALSE)
  }

  nv <- as.integer(graph_num_vertices(graph))
  ig <- igraph::make_empty_graph(n = nv, directed = isTRUE(graph$directed))

  edges <- graph_edges(graph)
  if (nrow(edges) > 0L) {
    ig <- igraph::add_edges(ig, as.vector(t(edges)))
  }

  if (include_properties) {
    for (nm in graph_graph_properties(graph)) {
      ig <- igraph::set_graph_attr(ig, nm, graph_graph_property(graph, nm))
    }
    for (nm in graph_vertex_properties(graph)) {
      val <- graph_vertex_property(graph, nm)
      if (is.matrix(val)) {
        val <- lapply(seq_len(nrow(val)), function(i) unname(val[i, ]))
      }
      ig <- igraph::set_vertex_attr(ig, nm, value = val)
    }
    for (nm in graph_edge_properties(graph)) {
      val <- graph_edge_property(graph, nm)
      if (is.matrix(val)) {
        val <- lapply(seq_len(nrow(val)), function(i) unname(val[i, ]))
      }
      ig <- igraph::set_edge_attr(ig, nm, value = val)
    }
  }

  ig
}

#' @noRd
#' @export
from_igraph <- function(graph, include_properties = TRUE) {
  require_igraph("from_igraph")
  if (!inherits(graph, "igraph")) {
    stop("graph must be an igraph object", call. = FALSE)
  }
  if (!is.logical(include_properties) || length(include_properties) != 1L ||
      is.na(include_properties)) {
    stop("include_properties must be TRUE or FALSE", call. = FALSE)
  }

  nv <- as.integer(igraph::vcount(graph))
  directed <- isTRUE(igraph::is_directed(graph))
  g <- graph_create(nv, directed = directed)

  edge_props <- NULL
  if (include_properties) {
    enms <- igraph::edge_attr_names(graph)
    if (length(enms)) {
      edge_props <- stats::setNames(
        lapply(enms, function(nm) igraph::edge_attr(graph, nm)),
        enms
      )
    }
  }

  el <- igraph::as_edgelist(graph, names = FALSE)
  if (nrow(el) > 0L) {
    graph_add_edges(g, el, properties = edge_props)
  } else if (!is.null(edge_props)) {
    for (nm in names(edge_props)) {
      graph_set_edge_property(g, nm, edge_props[[nm]])
    }
  }

  if (include_properties) {
    for (nm in igraph::graph_attr_names(graph)) {
      graph_set_graph_property(g, nm, igraph::graph_attr(graph, nm))
    }
    for (nm in igraph::vertex_attr_names(graph)) {
      graph_set_vertex_property(g, nm, igraph::vertex_attr(graph, nm))
    }
  }

  g
}

infer_graph_format <- function(file) {
  lf <- tolower(file)
  if (grepl("\\.gt(\\.(gz|bz2|xz|zst|zstd))?$", lf)) return("gt")
  if (grepl("\\.(graphml|xml)(\\.gz)?$", lf)) return("graphml")
  if (grepl("\\.gml$", lf)) return("gml")
  if (grepl("\\.(edgelist|ncol|lgl)$", lf)) return("edgelist")
  stop(sprintf("Cannot infer graph format from filename '%s'; specify fmt explicitly", file),
       call. = FALSE)
}

#' @noRd
#' @export
graph_load <- function(file, fmt = c("auto", "gt", "graphml", "gml", "edgelist"), ...) {
  fmt <- match.arg(fmt)
  if (fmt == "auto") {
    fmt <- infer_graph_format(file)
  }
  if (fmt == "gt") {
    return(read_gt(file))
  }
  require_igraph("graph_load")
  ig <- igraph::read_graph(path.expand(file), format = fmt, ...)
  from_igraph(ig)
}

#' @noRd
#' @export
graph_save <- function(graph, file, fmt = c("auto", "gt", "graphml", "gml", "edgelist"), ...) {
  graph_pointer(graph)
  fmt <- match.arg(fmt)
  if (fmt == "auto") {
    fmt <- infer_graph_format(file)
  }
  if (fmt == "gt") {
    return(write_gt(graph, file, ...))
  }
  require_igraph("graph_save")
  ig <- as_igraph(graph)
  igraph::write_graph(ig, path.expand(file), format = fmt, ...)
  invisible(file)
}
