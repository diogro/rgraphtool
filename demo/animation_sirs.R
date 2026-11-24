# Cookbook demo: SIRS epidemic animation (and fear-state variant)
# R-native counterpart of graph-tool's animation_sirs.py and animation_zombies.py.
# Runs fully offline; exports interactive HTML5 Canvas animations and sample PNG
# frames into tempdir().

library(rgraphtool)

set.seed(42)

# Build a 2D spatial proximity graph (k-nearest-neighbor mesh on random 2D points),
# analogous to the Delaunay triangulation in graph-tool's animation_sirs.py.
n <- 120L
pts <- cbind(x = stats::runif(n), y = stats::runif(n))
dmat <- as.matrix(stats::dist(pts))
diag(dmat) <- Inf
k_nn <- 5L
edge_pairs <- do.call(rbind, lapply(seq_len(n), function(i) {
  nbrs <- order(dmat[i, ])[seq_len(k_nn)]
  cbind(pmin(i, nbrs), pmax(i, nbrs))
}))
edge_pairs <- unique(edge_pairs)

g <- graph_create(n, directed = FALSE)
graph_add_edges(g, edge_pairs)
pos <- sfdp_layout(g, pos = pts, K = 0.25, max_iter = 60L)

# Seed a few initially infected vertices so waves propagate immediately
init_state <- integer(n)
init_state[c(1L, 40L, 80L)] <- 1L

# 1. Standard SIRS animation with recovered vertex hiding and newly-infected halos
html_sirs <- file.path(tempdir(), "rgraphtool_demo_sirs.html")
anim_sirs <- animate_sirs(
  g,
  pos = pos,
  n_steps = 80L,
  x = 0.002,
  r = 0.12,
  s = 0.02,
  init_state = init_state,
  hide_recovered = TRUE,
  fear_state = FALSE,
  mode = "html",
  output = html_sirs,
  fps = 12,
  open = FALSE,
  main = "SIRS Epidemic Waves on a Spatial Network"
)
print(anim_sirs)

# 2. Fear-state ("zombie") variant highlighting susceptible neighbors of infected nodes
html_fear <- file.path(tempdir(), "rgraphtool_demo_sirs_fear.html")
set.seed(42)
anim_fear <- animate_sirs(
  g,
  pos = pos,
  n_steps = 80L,
  x = 0.002,
  r = 0.12,
  s = 0.02,
  init_state = init_state,
  hide_recovered = FALSE,
  fear_state = TRUE,
  mode = "html",
  output = html_fear,
  fps = 12,
  open = FALSE,
  main = "SIRS Epidemic with Exposed/Fear State"
)
print(anim_fear)

# 3. Offscreen PNG frame sequence (short sample ready for ffmpeg or gifski)
frames_dir <- file.path(tempdir(), "rgraphtool_demo_sirs_frames")
set.seed(42)
anim_png <- animate_sirs(
  g,
  pos = pos,
  n_steps = 6L,
  x = 0.002,
  r = 0.12,
  s = 0.02,
  init_state = init_state,
  hide_recovered = TRUE,
  mode = "frames",
  output = frames_dir,
  width = 480L,
  height = 360L
)
cat(sprintf("Wrote %d PNG frames to: %s\n", length(anim_png$files), frames_dir))
cat(sprintf("Open interactive HTML5 animations in any browser:\n  %s\n  %s\n",
            html_sirs, html_fear))
