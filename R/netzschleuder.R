read_netzschleuder_csv <- function(archive, member) {
  connection <- unz(archive, member)
  header <- readLines(connection, n = 1L, warn = FALSE)
  close(connection)
  if (!length(header) || !startsWith(header, "#")) {
    stop("Invalid Netzschleuder archive: ", member, " has no header line", call. = FALSE)
  }
  columns <- trimws(strsplit(sub("^#\\s*", "", header), ",", fixed = TRUE)[[1]])
  utils::read.csv(unz(archive, member), header = FALSE, skip = 1L, col.names = columns,
                  stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
}

parse_positions <- function(x) {
  if (!is.character(x) || any(!grepl("^array\\(\\[.*\\]\\)$", x))) {
    return(NULL)
  }
  coordinates <- strsplit(substring(x, 8L, nchar(x) - 2L), ",", fixed = TRUE)
  if (any(lengths(coordinates) != 2L)) return(NULL)
  positions <- suppressWarnings(matrix(as.numeric(unlist(coordinates)), ncol = 2L,
                                       byrow = TRUE))
  if (any(!is.finite(positions))) return(NULL)
  colnames(positions) <- c("x", "y")
  positions
}

netzschleuder_network <- function(name, network = name, file = NULL, directed = NULL) {
  valid_name <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[A-Za-z0-9_.-]+$", x)
  }
  if (!valid_name(name) || !valid_name(network)) {
    stop("name and network must be Netzschleuder identifiers", call. = FALSE)
  }
  if (is.null(file)) {
    file <- tempfile(fileext = ".zip")
    on.exit(unlink(file), add = TRUE)
    utils::download.file(
      sprintf("https://networks.skewed.de/net/%s/files/%s.csv.zip", name, network),
      file, mode = "wb", quiet = TRUE
    )
  } else if (!is.character(file) || length(file) != 1L || is.na(file) ||
             !file.exists(file) || dir.exists(file)) {
    stop("file must be the path to a Netzschleuder CSV archive", call. = FALSE)
  }
  members <- utils::unzip(file, list = TRUE)$Name
  if (!all(c("nodes.csv", "edges.csv") %in% members)) {
    stop("Invalid Netzschleuder archive: expected nodes.csv and edges.csv", call. = FALSE)
  }
  nodes <- read_netzschleuder_csv(file, "nodes.csv")
  edges <- read_netzschleuder_csv(file, "edges.csv")
  properties <- list()
  if ("gprops.csv" %in% members) {
    gprops <- read_netzschleuder_csv(file, "gprops.csv")
    properties <- stats::setNames(as.list(gprops[[2]]), gprops[[1]])
  }
  if (is.null(directed)) {
    directed <- !is.null(properties$directed) &&
      isTRUE(as.numeric(properties$directed) == 1)
  }
  n <- nrow(nodes)
  if (!"index" %in% names(nodes) || !is.numeric(nodes$index) ||
      !identical(sort(as.double(nodes$index)), as.double(seq_len(n) - 1L)) ||
      ncol(edges) < 2L || !is.numeric(edges[[1]]) || !is.numeric(edges[[2]])) {
    stop("Invalid Netzschleuder archive: malformed node indices or edge list", call. = FALSE)
  }
  nodes <- nodes[order(nodes$index), , drop = FALSE]
  nodes$index <- as.integer(nodes$index) + 1L
  positions <- NULL
  if ("_pos" %in% names(nodes)) {
    positions <- parse_positions(nodes[["_pos"]])
    nodes[["_pos"]] <- NULL
  }
  names(edges)[1:2] <- c("source", "target")
  edges$source <- edges$source + 1
  edges$target <- edges$target + 1
  graph <- graph_create(n, directed = directed)
  extra_cols <- setdiff(names(edges), c("source", "target"))
  num_props <- list()
  edge_levels <- list()
  for (col in extra_cols) {
    if (is.numeric(edges[[col]])) {
      num_props[[col]] <- as.numeric(edges[[col]])
    } else if (is.character(edges[[col]]) || is.factor(edges[[col]]) ||
               is.logical(edges[[col]])) {
      vals <- if (is.character(edges[[col]])) trimws(edges[[col]]) else edges[[col]]
      f <- factor(vals)
      num_props[[col]] <- as.numeric(as.integer(f))
      edge_levels[[col]] <- levels(f)
    }
  }
  graph_add_edges(graph, as.matrix(edges[, c("source", "target")]),
                  properties = if (length(num_props)) num_props else NULL)
  rownames(nodes) <- NULL
  res <- list(graph = graph, nodes = nodes, edges = edges, positions = positions,
              properties = properties)
  if (length(edge_levels)) res$edge_levels <- edge_levels
  res
}
