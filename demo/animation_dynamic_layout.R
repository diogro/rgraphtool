# Cookbook demo: Dynamic force-directed layout with spatial-segregation rewiring
# R-native counterpart of graph-tool's animation_dancing.py.
# Runs fully offline; exports an interactive HTML5 Canvas animation and sample
# PNG frames into tempdir().

library(rgraphtool)

set.seed(2026)

# Create a connected random graph (ring backbone + random chords) with 80 vertices
n <- 80L
ring <- cbind(seq_len(n), c(seq_len(n - 1L) + 1L, 1L))
chords <- matrix(0L, nrow = 120L, ncol = 2L)
for (i in seq_len(nrow(chords))) {
  u <- sample.int(n, 1L)
  v <- sample.int(n, 1L)
  while (v == u) v <- sample.int(n, 1L)
  chords[i, ] <- c(min(u, v), max(u, v))
}
edges <- unique(rbind(ring, chords))

g <- graph_create(n, directed = FALSE)
graph_add_edges(g, edges)

# Initial SFDP force-directed layout
pos0 <- sfdp_layout(g, K = 0.5, max_iter = 80L)

# Run incremental SFDP + distance-shortening rewiring (spatial segregation)
html_dyn <- file.path(tempdir(), "rgraphtool_demo_dynamic_layout.html")
dyn_res <- animate_dynamic_layout(
  g,
  pos = pos0,
  n_steps = 80L,
  K = 0.5,
  step = 0.005,
  n_rewire = 25L,
  mode = "html",
  output = html_dyn,
  fps = 15,
  open = FALSE,
  main = "Dynamic Layout & Spatial Segregation"
)
print(dyn_res)

pct_drop <- 100 * (1 - dyn_res$mean_edge_length[length(dyn_res$mean_edge_length)] /
                     dyn_res$mean_edge_length[1L])
cat(sprintf("Mean edge length contracted by %.1f%% over %d steps (%d accepted rewirings).\n",
            pct_drop, dyn_res$params$n_steps, sum(dyn_res$rewired)))

# Also export a short sequence of offscreen PNG frames
frames_dir <- file.path(tempdir(), "rgraphtool_demo_dynlayout_frames")
set.seed(2026)
dyn_png <- animate_dynamic_layout(
  g,
  pos = pos0,
  n_steps = 5L,
  K = 0.5,
  step = 0.005,
  n_rewire = 25L,
  mode = "frames",
  output = frames_dir,
  width = 480L,
  height = 360L
)
cat(sprintf("Wrote %d PNG frames to: %s\n", length(dyn_png$files), frames_dir))
cat(sprintf("Open interactive HTML5 animation in any browser:\n  %s\n", html_dyn))
