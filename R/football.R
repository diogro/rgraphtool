football_network <- function(file = NULL) {
  if (is.null(file)) {
    file <- tempfile(fileext = ".zip")
    on.exit(unlink(file), add = TRUE)
    utils::download.file(
      "https://networks.skewed.de/net/football/files/football.csv.zip",
      file, mode = "wb", quiet = TRUE
    )
  } else if (!is.character(file) || length(file) != 1L || is.na(file) ||
             !file.exists(file) || dir.exists(file)) {
    stop("file must be the path to a football CSV archive", call. = FALSE)
  }

  nodes <- utils::read.csv(unz(file, "nodes.csv"), header = FALSE,
                           comment.char = "#", stringsAsFactors = FALSE,
                           col.names = c("index", "label", "conference", "position"))
  edges <- utils::read.csv(unz(file, "edges.csv"), header = FALSE,
                           comment.char = "#", col.names = c("source", "target"))
  if (nrow(nodes) != 115L || nrow(edges) != 613L ||
      anyNA(nodes) || anyNA(edges) || !is.numeric(nodes$index) ||
      anyDuplicated(nodes$index) ||
      !identical(sort(as.double(nodes$index)), as.double(seq_len(115L) - 1L))) {
    stop("Invalid football archive: expected 115 indexed nodes and 613 edges",
         call. = FALSE)
  }
  nodes <- nodes[order(nodes$index), , drop = FALSE]
  if (any(!grepl("^array\\(\\[.*\\]\\)$", nodes$position))) {
    stop("Invalid football archive: malformed node positions", call. = FALSE)
  }
  coordinates <- strsplit(substring(nodes$position, 8L, nchar(nodes$position) - 2L),
                          ",", fixed = TRUE)
  if (any(lengths(coordinates) != 2L)) {
    stop("Invalid football archive: malformed node positions", call. = FALSE)
  }
  positions <- suppressWarnings(matrix(as.numeric(unlist(coordinates)), ncol = 2L, byrow = TRUE))
  if (any(!is.finite(positions))) {
    stop("Invalid football archive: non-numeric node positions", call. = FALSE)
  }
  colnames(positions) <- c("x", "y")
  endpoints <- as.matrix(edges)
  graph <- graph_create(nrow(nodes), directed = FALSE)
  graph_add_edges(graph, endpoints + 1)
  nodes$index <- nodes$index + 1L
  nodes$position <- NULL
  rownames(nodes) <- NULL
  list(graph = graph, nodes = nodes, positions = positions)
}
