# Cookbook demo: Interactive BFS shortest-path tree visualization
# R-native counterpart of graph-tool's interactive_bst.py.
# Runs fully offline on a synthetic modular collaboration network; generates a
# zero-dependency interactive HTML5 Canvas explorer and a multi-root PNG snapshot.

library(rgraphtool)

set.seed(7)

# Build a connected modular collaboration graph (6 research groups x 15 nodes = 90)
n_groups <- 6L
group_size <- 15L
n <- n_groups * group_size
groups <- rep(seq_len(n_groups), each = group_size)

edge_list <- list()
for (k in seq_len(n_groups)) {
  nodes <- which(groups == k)
  # Local ring inside each group for guaranteed connectivity
  edge_list[[length(edge_list) + 1L]] <- cbind(nodes, c(nodes[-1L], nodes[1L]))
  # Random intra-group collaboration edges
  for (i in seq_len(group_size - 1L)) {
    for (j in (i + 1L):group_size) {
      if (stats::runif(1L) < 0.22) {
        edge_list[[length(edge_list) + 1L]] <- cbind(nodes[i], nodes[j])
      }
    }
  }
}
# Sparse inter-group bridges connecting neighboring and distant groups
for (k in seq_len(n_groups)) {
  k2 <- if (k < n_groups) k + 1L else 1L
  u <- sample(which(groups == k), 2L)
  v <- sample(which(groups == k2), 2L)
  edge_list[[length(edge_list) + 1L]] <- cbind(pmin(u, v), pmax(u, v))
}
edges <- unique(do.call(rbind, edge_list))

g <- graph_create(n, directed = FALSE)
graph_add_edges(g, edges)
pos <- sfdp_layout(g, groups = groups, mu = 0.35, max_iter = 100L)

# 1. Interactive HTML5 Canvas explorer (hover or click any vertex in a browser)
html_bfs <- file.path(tempdir(), "rgraphtool_demo_interactive_bfs.html")
ib <- interactive_bfs(
  g,
  pos = pos,
  source = 1L,
  max_dist = 3L,
  mode = "html",
  output = html_bfs,
  open = FALSE,
  main = "Interactive BFS Tree Explorer (Modular Network)"
)
print(ib)

# 2. Multi-root static PNG snapshot comparing depth-3 BFS trees from 4 vertices
snap_png <- file.path(tempdir(), "rgraphtool_demo_bfs_snapshot.png")
ib_snap <- interactive_bfs(
  g,
  pos = pos,
  source = c(1L, 20L, 48L, 76L),
  max_dist = 3L,
  mode = "snapshot",
  output = snap_png,
  width = 800L,
  height = 640L,
  main = "BFS Tree"
)
cat(sprintf("Wrote multi-root BFS snapshot PNG: %s\n", snap_png))
cat(sprintf("Open interactive HTML5 BFS explorer in any browser:\n  %s\n", html_bfs))
cat("In an interactive R session, call interactive_bfs(g, pos, mode = 'locator') for click-to-inspect in R graphics.\n")
