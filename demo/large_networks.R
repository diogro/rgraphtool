# Large-network stochastic block model inference: strictly agglomerative
# search, hierarchy-shape initialization (simple_init), warm-starting, and
# multi-threaded execution (graph-tool cookbook, "Large networks").
# Generates a reproducible 10,000-vertex, ~50,000-edge modular network so the
# demo runs completely offline in a few seconds.
library(rgraphtool)

set.seed(42)
n <- 10000L
B_true <- 50L
k_avg <- 10
sz <- n %/% B_true
m_in <- as.integer(round(n * k_avg * 0.42))
m_out <- as.integer(round(n * k_avg * 0.08))
g_id <- sample.int(B_true, m_in, replace = TRUE)
u_in <- (g_id - 1L) * sz + sample.int(sz, m_in, replace = TRUE)
v_in <- (g_id - 1L) * sz + sample.int(sz, m_in, replace = TRUE)
keep <- u_in != v_in & u_in <= n & v_in <= n
u_out <- sample.int(n, m_out, replace = TRUE)
v_out <- sample.int(n, m_out, replace = TRUE)
keep2 <- u_out != v_out
edges <- rbind(cbind(u_in[keep], v_in[keep]), cbind(u_out[keep2], v_out[keep2]))
g <- graph_add_edges(graph_create(n, directed = FALSE), edges)
truth <- rep(seq_len(B_true), each = sz)[seq_len(n)]

nthr <- min(4L, rgraphtool_max_threads())
cat(sprintf(
  "Synthetic large network: N = %d vertices, E = %d edges, true B = %d (using up to %d threads)\n",
  graph_num_vertices(g), graph_num_edges(g), B_true, nthr
))

# Warm up the partition-count cache once before timing
invisible(block_description_length(block_state(g)))

# 1. Strictly agglomerative flat SBM (bisection = FALSE, r = 0.5): 1 thread vs multi-thread
set.seed(1)
t0 <- proc.time()[3]
st_1t <- minimize_blockmodel_dl(
  g, nthreads = 1L, control = list(bisection = FALSE, r = 0.5)
)
t_flat_1t <- proc.time()[3] - t0

set.seed(1)
t0 <- proc.time()[3]
st_mt <- minimize_blockmodel_dl(
  g, nthreads = nthr, control = list(bisection = FALSE, r = 0.5)
)
t_flat_mt <- proc.time()[3] - t0

cat(sprintf(
  "Flat agglom (1 thread):   %.3f s | B = %d | DL = %.1f | overlap = %.4f\n",
  t_flat_1t, length(block_sizes(st_1t)), block_description_length(st_1t),
  partition_overlap(block_membership(st_1t), truth, norm = TRUE)
))
cat(sprintf(
  "Flat agglom (%d threads):  %.3f s | B = %d | DL = %.1f | overlap = %.4f (%.2fx speedup)\n",
  nthr, t_flat_mt, length(block_sizes(st_mt)), block_description_length(st_mt),
  partition_overlap(block_membership(st_mt), truth, norm = TRUE),
  t_flat_1t / max(t_flat_mt, 1e-6)
))

# 2. Nested SBM with hierarchy-shape initialization (simple_init = TRUE)
set.seed(1)
t0 <- proc.time()[3]
nst_1t <- minimize_nested_blockmodel_dl(
  g, simple_init = TRUE, B_min_base = 200L, nthreads = 1L, control = list(r = 0.5)
)
t_nest_1t <- proc.time()[3] - t0

set.seed(1)
t0 <- proc.time()[3]
nst_mt <- minimize_nested_blockmodel_dl(
  g, simple_init = TRUE, B_min_base = 200L, nthreads = nthr, control = list(r = 0.5)
)
t_nest_mt <- proc.time()[3] - t0

cat(sprintf(
  "Nested simple_init (1t):  %.3f s | B0 = %d | DL = %.1f | overlap = %.4f\n",
  t_nest_1t, max(nested_levels(nst_1t)[[1L]]), nested_description_length(nst_1t),
  partition_overlap(nested_levels(nst_1t)[[1L]], truth, norm = TRUE)
))
cat(sprintf(
  "Nested simple_init (%dt):  %.3f s | B0 = %d | DL = %.1f | overlap = %.4f (%.2fx speedup)\n",
  nthr, t_nest_mt, max(nested_levels(nst_mt)[[1L]]), nested_description_length(nst_mt),
  partition_overlap(nested_levels(nst_mt)[[1L]], truth, norm = TRUE),
  t_nest_1t / max(t_nest_mt, 1e-6)
))
cat("\nRecovered nested hierarchy summary:\n")
print(nested_summary(nst_mt))

# 3. Warm-starting refinement from an existing hierarchy
t0 <- proc.time()[3]
nst_ref <- minimize_nested_blockmodel_dl(g, init_levels = nst_mt, nthreads = nthr)
cat(sprintf("Warm-start nested refine: %.3f s | DL = %.1f\n",
            proc.time()[3] - t0, nested_description_length(nst_ref)))

# 4. Visualization: runtime comparison and recovered community sizes
old_par <- graphics::par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
on.exit(graphics::par(old_par), add = TRUE)

times <- c(t_flat_1t, t_flat_mt, t_nest_1t, t_nest_mt)
labels <- c("Flat\n(1 thr)", sprintf("Flat\n(%d thr)", nthr),
            "Nested\n(1 thr)", sprintf("Nested\n(%d thr)", nthr))
graphics::barplot(
  times, names.arg = labels, ylab = "Wall time (s)",
  col = c("grey60", "steelblue", "grey60", "darkorange"),
  main = sprintf("Large-network scaling (N=%d, E=%d)", n, graph_num_edges(g))
)

graphics::hist(
  block_sizes(st_mt), breaks = 15, col = "steelblue", border = "white",
  xlab = "Community size (vertices)",
  main = sprintf("Recovered community sizes (B = %d)", length(block_sizes(st_mt)))
)
