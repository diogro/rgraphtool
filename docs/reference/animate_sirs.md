# SIRS epidemic and dynamic force-directed layout animations

R-native equivalents of graph-tool's SIRS epidemic animation
(`animation_sirs.py`, `animation_zombies.py`) and dynamic force-directed
spatial-segregation animation (`animation_dancing.py`). Simulation loops
execute natively in C++ and respect
[`set.seed`](https://rdrr.io/r/base/Random.html).

## Usage

``` r
animate_sirs(graph, pos = NULL, n_steps = 120L, x = 0.001, r = 0.1, s = 0.01,
  init_state = NULL, directed = FALSE, hide_recovered = TRUE,
  fear_state = FALSE, mode = c("html", "frames", "device", "none"),
  output = NULL, width = 760L, height = 560L, fps = 12, delay = 0.05,
  open = FALSE, main = "SIRS Epidemic Animation")

animate_dynamic_layout(graph, pos = NULL, n_steps = 80L, K = 0.5, C = 0.2,
  p = 2.0, step = 0.005, n_rewire = 10L,
  mode = c("html", "frames", "device", "none"), output = NULL, width = 760L,
  height = 560L, fps = 15, delay = 0.05, open = FALSE,
  main = "Dynamic Layout & Spatial Segregation")
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- pos:

  Optional \\N \times 2\\ initial vertex layout matrix. Computed with
  [`sfdp_layout`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md)
  if `NULL`.

- n_steps:

  Number of simulation steps to run (produces `n_steps + 1` frames
  including step 0).

- x:

  Spontaneous infection probability (`S -> I`) per vertex per step.

- r:

  Recovery probability (`I -> R`) per vertex per step.

- s:

  Resusceptibility probability (`R -> S`) per vertex per step.

- init_state:

  Optional integer vector of length \\N\\ with values in `{0, 1, 2}`
  (`0 = S`, `1 = I`, `2 = R`). Defaults to all susceptible.

- directed:

  Logical; if `TRUE` and `graph$directed` is `TRUE`, infection spreads
  along directed edges `u -> v`; otherwise undirected neighbors are
  used.

- hide_recovered:

  Logical; if `TRUE`, recovered vertices and their incident edges are
  hidden from the animation (matching graph-tool's `animation_sirs.py`).

- fear_state:

  Logical; if `TRUE`, susceptible vertices adjacent to at least one
  infected neighbor are highlighted in amber (matching graph-tool's
  `animation_zombies.py`).

- K, C, p:

  SFDP spring-electrical layout parameters for
  `animate_dynamic_layout()`.

- step:

  Per-frame incremental SFDP displacement step size.

- n_rewire:

  Number of distance-shortening edge rewiring attempts per frame
  (preserving edge count and avoiding isolated vertices, self-loops, and
  parallel edges).

- mode:

  Output mode: `"html"` writes a self-contained zero-dependency
  interactive HTML5 Canvas viewer; `"frames"` writes numbered PNG files
  to a directory; `"device"` animates on the active R graphics device;
  `"none"` returns the trajectory without rendering.

- output:

  Optional destination HTML file path (`mode = "html"`) or directory
  path (`mode = "frames"`). Defaults to a temporary path.

- width, height:

  Canvas or PNG frame dimensions in pixels.

- fps:

  Initial playback rate (frames per second) for `mode = "html"`.

- delay:

  Pause in seconds between frames when `mode = "device"`.

- open:

  Logical; if `TRUE` in an interactive R session with `mode = "html"`,
  opens the generated HTML file via
  [`browseURL`](https://rdrr.io/r/utils/browseURL.html).

- main:

  Animation title.

## Details

**Platform and display requirements:**

- `mode = "html"` requires no display server, GTK+, Cairo, or external R
  packages; the generated standalone HTML5 file runs offline in any
  modern web browser or RStudio Viewer with interactive play/pause,
  frame scrubbing, speed control, and live summary charts.

- `mode = "frames"` requires R's built-in
  [`png`](https://rdrr.io/r/grDevices/png.html) device (available on
  standard macOS, Windows, and Cairo/X11-enabled Linux builds) and
  produces numbered PNG frames suitable for external encoders such as
  `ffmpeg`.

- `mode = "device"` renders directly to the active interactive R
  graphics device (Quartz, X11, Windows, or RStudioGD).

## Value

`animate_sirs()` returns an object of class
`"rgraphtool_sirs_animation"` containing `pos`, `states`,
`newly_infected`, `fear`, `counts`, `params`, `mode`, `output`, and
`files`.

`animate_dynamic_layout()` returns an object of class
`"rgraphtool_dynamic_layout"` containing `positions`, `edges`,
`mean_edge_length`, `rewired`, `params`, `mode`, `output`, and `files`.

## Examples

``` r
g <- graph_create(20L, directed = FALSE)
graph_add_edges(g, cbind(1:20, c(2:20, 1L)))
set.seed(1)
init <- integer(20L); init[1L] <- 1L
sirs <- animate_sirs(g, n_steps = 10L, init_state = init, mode = "none")
print(sirs)
#> rgraphtool SIRS epidemic animation
#>   vertices: 20 | frames: 11 (steps = 10)
#>   rates: x = 0.001, r = 0.1, s = 0.01 | hide_recovered = TRUE | fear_state = FALSE
#>   peak infected: 7 | total new infections: 10
#>   mode: none
dyn <- animate_dynamic_layout(g, n_steps = 10L, n_rewire = 5L, mode = "none")
print(dyn)
#> rgraphtool dynamic layout animation
#>   vertices: 20 | edges: 20 | frames: 11
#>   mean edge length: 0.4440 -> 0.4440 | total rewirings: 0
#>   mode: none
```
