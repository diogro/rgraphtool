# Cookbook demo: Disentangling homophily / community structure and triadic
# closure using LatentClosureBlockState (Peixoto 2022 PRX).

library(rgraphtool)
set.seed(2026)

cat("=== Disentangling substrate edges vs. triadic closure edges ===\n")

# Build a sparse 24-vertex network:
#   - Substrate: two 12-cycles connected by 3 bridge edges (triangle-free!)
#   - Triadic closures: 12 chords closing open 2-paths (i, i+2) through i+1
N <- 24L
cyc1 <- rbind(cbind(1:11, 2:12), c(12L, 1L))
cyc2 <- rbind(cbind(13:23, 14:24), c(24L, 13L))
bridges <- matrix(c(1L, 13L, 5L, 17L, 9L, 21L), ncol = 2L, byrow = TRUE)
sub_edges <- rbind(cyc1, cyc2, bridges)

cl1 <- cbind(seq(1L, 11L, by = 2L), (seq(3L, 13L, by = 2L) - 1L) %% 12L + 1L)
cl2 <- cbind(seq(13L, 23L, by = 2L), 12L + ((seq(3L, 13L, by = 2L) - 1L) %% 12L + 1L))
cl_edges <- rbind(cl1, cl2)

g <- graph_create(N, directed = FALSE)
graph_add_edges(g, rbind(sub_edges, cl_edges))

cat(sprintf(
  "Graph: %d vertices, %d edges (%d substrate + %d true triadic closures)\n",
  graph_num_vertices(g), graph_num_edges(g), nrow(sub_edges), nrow(cl_edges)
))

st <- latent_closure_block_state(g, L = 1L, nested = TRUE)
cat(sprintf("Initial DL (all edges in substrate): %.3f\n", closure_description_length(st)))

res <- infer_triadic_closure(st, nsamples = 250L, burnin = 200L, thin = 2L)
print(res)

ek <- function(u, v) paste0(pmin(u, v), ":", pmax(u, v))
is_closure <- ek(res$edges$from, res$edges$to) %in% ek(cl_edges[, 1L], cl_edges[, 2L])

cat(sprintf(
  "Posterior mean closure_prob: true closure edges = %.3f vs. substrate edges = %.3f\n",
  mean(res$edges$closure_prob[is_closure]),
  mean(res$edges$closure_prob[!is_closure])
))

ev <- evaluate_edge_prediction(res$edges$closure_prob, as.integer(is_closure))
cat(sprintf("Closure vs. substrate separation ROC-AUC = %.4f\n\n", ev$roc_auc))

vs <- res$vertex_stats$gen_1
cat("Top 8 closing ego vertices (t = closed triads, m = open triads, c = fraction):\n")
print(utils::head(vs[order(-vs$t), ], 8L), row.names = FALSE)
