test_that("edge property maps preserve edge alignment across unsorted edges and updates", {
  g <- graph_create(5, directed = FALSE)
  # Add edges in non-sorted source order so boost::edges order differs from insertion order
  e1 <- rbind(c(4L, 5L), c(1L, 2L), c(3L, 5L), c(1L, 4L))
  g <- graph_add_edges(g, e1, properties = list(w = c(40, 10, 30, 20), tag = c(4, 1, 3, 2)))
  expect_setequal(graph_edge_properties(g), c("w", "tag"))

  ed <- graph_edges(g)
  w <- graph_edge_property(g, "w")
  expect_equal(length(w), nrow(ed))
  # Match each edge back to its expected weight regardless of boost::edges order
  key <- function(m) paste(pmin(m[, 1], m[, 2]), pmax(m[, 1], m[, 2]), sep = "-")
  expected <- c("4-5" = 40, "1-2" = 10, "3-5" = 30, "1-4" = 20)
  expect_equal(unname(w), unname(expected[key(ed)]))

  # Adding more edges without properties pads existing properties with NA; with properties fills them
  g <- graph_add_edges(g, rbind(c(2L, 3L)), properties = list(w = 50))
  ed2 <- graph_edges(g)
  w2 <- graph_edge_property(g, "w")
  tag2 <- graph_edge_property(g, "tag")
  expected2 <- c(expected, "2-3" = 50)
  expect_equal(unname(w2), unname(expected2[key(ed2)]))
  expect_true(is.na(tag2[key(ed2) == "2-3"]))

  # Setting a new property in graph_edges(g) order
  g <- graph_set_edge_property(g, "doubled", 2 * w2)
  expect_equal(graph_edge_property(g, "doubled"), 2 * w2)

  expect_error(graph_edge_property(g, "missing_prop"), "does not exist")
  expect_error(graph_set_edge_property(g, "bad_len", c(1, 2)), "length")
})

test_that("netzschleuder_network attaches numeric edge columns as edge properties", {
  td <- tempfile("ns_weighted_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  edges_csv <- file.path(td, "edges.csv")
  nodes_csv <- file.path(td, "nodes.csv")
  writeLines(c("# source, target, value, weight",
               "0, 1, 5, 1.5",
               "1, 2, 2, 0.4",
               "2, 0, 8, 3.2"), edges_csv)
  writeLines(c("# index, name",
               "0, A",
               "1, B",
               "2, C"), nodes_csv)
  zfile <- file.path(td, "toy_weighted.csv.zip")
  oldwd <- setwd(td)
  on.exit(setwd(oldwd), add = TRUE)
  utils::zip(zfile, files = c("edges.csv", "nodes.csv"), flags = "-q")
  setwd(oldwd)

  net <- netzschleuder_network("toy_weighted", file = zfile)
  expect_setequal(graph_edge_properties(net$graph), c("value", "weight"))
  expect_equal(sort(graph_edge_property(net$graph, "value")), c(2, 5, 8))
  expect_equal(sort(graph_edge_property(net$graph, "weight")), c(0.4, 1.5, 3.2))
})

test_that("transform_weights computes exact transformed values and log-Jacobians", {
  x <- c(1.5, 3.0, 6.2, 12.4)

  y_id <- transform_weights(x, "identity")
  expect_equal(as.numeric(y_id), x)
  expect_equal(attr(y_id, "log_jacobian"), 0)

  y_log <- transform_weights(x, "log")
  expect_equal(as.numeric(y_log), log(x))
  expect_equal(attr(y_log, "log_jacobian"), -sum(log(x)), tolerance = 1e-12)

  y_l1p <- transform_weights(x, "log1p")
  expect_equal(as.numeric(y_l1p), log1p(x))
  expect_equal(attr(y_l1p, "log_jacobian"), -sum(log1p(x)), tolerance = 1e-12)

  y_sqrt <- transform_weights(x, "sqrt")
  expect_equal(as.numeric(y_sqrt), sqrt(x))
  expect_equal(attr(y_sqrt, "log_jacobian"), sum(log(0.5 / sqrt(x))), tolerance = 1e-12)

  y_z <- transform_weights(x, "zscore")
  expect_equal(mean(y_z), 0, tolerance = 1e-12)
  expect_equal(stats::sd(y_z), 1, tolerance = 1e-12)
  expect_equal(attr(y_z, "log_jacobian"), -length(x) * log(stats::sd(x)), tolerance = 1e-12)

  expect_equal(as.numeric(transform_weights(c(-2, 0, 5), "sign")), c(-1, 0, 1))
  expect_true(is.na(attr(transform_weights(c(-2, 0, 5), "sign"), "log_jacobian")))
  expect_equal(as.numeric(transform_weights(c(-2, 0, 5), "abs")), c(2, 0, 5))
  expect_equal(as.numeric(transform_weights(c(10, 2, 5), "rank")), c(3, 1, 2))

  expect_error(transform_weights(c(0, 1, 2), "log"), "positive")
  expect_error(transform_weights(c(-1, 1, 2), "sqrt"), "positive")
})

test_that("flat and nested weighted SBM description length matches pinned graph-tool 3.9 values", {
  edges <- rbind(
    c(1L, 2L), c(1L, 2L), c(1L, 3L), c(1L, 1L), c(2L, 3L),
    c(3L, 4L), c(3L, 4L), c(4L, 5L), c(4L, 6L), c(5L, 6L), c(5L, 5L)
  )
  b_flat <- c(1L, 1L, 2L, 2L, 3L, 3L)
  bs_nest <- list(c(1L, 1L, 2L, 2L, 3L, 3L), c(1L, 1L, 2L), c(1L, 1L))

  w_exp <- c(1.2, 2.5, 0.7, 3.1, 1.5, 0.9, 2.1, 0.4, 1.8, 0.5, 2.2)
  w_norm <- c(-1.2, 2.5, 0.7, 3.1, -0.5, 0.9, 1.4, -2.1, 0.3, 1.1, -0.8)
  w_geom <- c(0, 2, 1, 4, 0, 3, 2, 1, 5, 1, 0)
  w_pois <- c(1, 3, 0, 5, 2, 1, 4, 2, 0, 3, 1)
  w_bin <- c(1, 3, 0, 4, 2, 1, 3, 2, 0, 4, 1)

  # Exact reference values from pinned graph-tool 3.9 (DC, distributed degree prior)
  ref <- list(
    list(dir = FALSE, type = "real-exponential",   mode = "microcanonical", w = w_exp,
         flat = 42.246464181403, nodl = 34.115263524035, nest = 47.593708181195),
    list(dir = FALSE, type = "real-exponential",   mode = "bayesian",       w = w_exp,
         flat = 45.853025273659, nodl = 45.853025273659, nest = 59.331469930818),
    list(dir = FALSE, type = "real-normal",        mode = "microcanonical", w = w_norm,
         flat = 39.484382256956, nodl = 34.900089839281, nest = 47.793092595705),
    list(dir = FALSE, type = "real-normal",        mode = "bayesian",       w = w_norm,
         flat = 49.790091657028, nodl = 49.790091657028, nest = 73.391480867489),
    list(dir = FALSE, type = "discrete-geometric", mode = "microcanonical", w = w_geom,
         flat = 44.850692737620, nodl = 35.761955187403, nest = 49.860220699482),
    list(dir = FALSE, type = "discrete-geometric", mode = "bayesian",       w = w_geom,
         flat = 49.662394956152, nodl = 49.662394956152, nest = 63.760660468232),
    list(dir = FALSE, type = "discrete-poisson",   mode = "microcanonical", w = w_pois,
         flat = 46.598368680310, nodl = 36.985902101491, nest = 51.877039363446),
    list(dir = FALSE, type = "discrete-poisson",   mode = "bayesian",       w = w_pois,
         flat = 54.283125516778, nodl = 54.283125516778, nest = 69.174262778733),
    list(dir = FALSE, type = "discrete-binomial",  mode = "microcanonical", w = w_bin,
         flat = 47.379022163628, nodl = 37.933609669472, nest = 52.494789375827),
    list(dir = FALSE, type = "discrete-binomial",  mode = "bayesian",       w = w_bin,
         flat = 49.287457336279, nodl = 49.287457336279, nest = 63.848637042633),
    list(dir = TRUE,  type = "real-exponential",   mode = "microcanonical", w = w_exp,
         flat = 47.533264990623, nodl = 39.402064333254, nest = 52.774013914531),
    list(dir = TRUE,  type = "real-exponential",   mode = "bayesian",       w = w_exp,
         flat = 51.139826082878, nodl = 51.139826082878, nest = 64.511775664155),
    list(dir = TRUE,  type = "real-normal",        mode = "microcanonical", w = w_norm,
         flat = 44.771183066175, nodl = 40.186890648500, nest = 52.973398329042),
    list(dir = TRUE,  type = "real-normal",        mode = "bayesian",       w = w_norm,
         flat = 55.076892466248, nodl = 55.076892466248, nest = 78.571786600826),
    list(dir = TRUE,  type = "discrete-geometric", mode = "microcanonical", w = w_geom,
         flat = 50.137493546839, nodl = 41.048755996622, nest = 55.040526432819),
    list(dir = TRUE,  type = "discrete-geometric", mode = "bayesian",       w = w_geom,
         flat = 54.949195765371, nodl = 54.949195765371, nest = 68.940966201568),
    list(dir = TRUE,  type = "discrete-poisson",   mode = "microcanonical", w = w_pois,
         flat = 51.885169489529, nodl = 42.272702910710, nest = 57.057345096782),
    list(dir = TRUE,  type = "discrete-poisson",   mode = "bayesian",       w = w_pois,
         flat = 59.569926325997, nodl = 59.569926325997, nest = 74.354568512070),
    list(dir = TRUE,  type = "discrete-binomial",  mode = "microcanonical", w = w_bin,
         flat = 52.665822972847, nodl = 43.220410478691, nest = 57.675095109164),
    list(dir = TRUE,  type = "discrete-binomial",  mode = "bayesian",       w = w_bin,
         flat = 54.574258145498, nodl = 54.574258145498, nest = 69.028942775970)
  )

  g_u <- graph_add_edges(graph_create(6, directed = FALSE), edges)
  g_d <- graph_add_edges(graph_create(6, directed = TRUE), edges)

  for (r in ref) {
    g <- if (r$dir) g_d else g_u
    st <- weighted_block_state(g, weights = r$w, types = r$type, params = r$mode,
                               blocks = b_flat, degree_corrected = TRUE)
    expect_equal(weighted_description_length(st), r$flat, tolerance = 1e-10,
                 info = paste(r$dir, r$type, r$mode, "flat"))
    expect_equal(weighted_description_length(st, rec_dl = FALSE), r$nodl, tolerance = 1e-10,
                 info = paste(r$dir, r$type, r$mode, "nodl"))

    nst <- nested_weighted_block_state(g, weights = r$w, types = r$type, params = r$mode,
                                       levels = bs_nest, degree_corrected = TRUE)
    expect_equal(nested_weighted_description_length(nst), r$nest, tolerance = 1e-10,
                 info = paste(r$dir, r$type, r$mode, "nested"))
  }

  # Simultaneous multi-covariate check (real-normal + discrete-poisson)
  w_mat <- cbind(w_norm, w_pois)
  st_mu <- weighted_block_state(g_u, weights = w_mat,
                                types = c("real-normal", "discrete-poisson"),
                                blocks = b_flat)
  expect_equal(weighted_description_length(st_mu), 58.776964239442, tolerance = 1e-10)
  nst_mu <- nested_weighted_block_state(g_u, weights = w_mat,
                                        types = c("real-normal", "discrete-poisson"),
                                        levels = bs_nest)
  expect_equal(nested_weighted_description_length(nst_mu), 68.124046786482, tolerance = 1e-10)

  st_md <- weighted_block_state(g_d, weights = w_mat,
                                types = c("real-normal", "discrete-poisson"),
                                blocks = b_flat)
  expect_equal(weighted_description_length(st_md), 64.063765048661, tolerance = 1e-10)
  nst_md <- nested_weighted_block_state(g_d, weights = w_mat,
                                        types = c("real-normal", "discrete-poisson"),
                                        levels = bs_nest)
  expect_equal(nested_weighted_description_length(nst_md), 73.304352519818, tolerance = 1e-10)
})

test_that("incremental weighted virtual_move/move deltas match full entropy recomputation", {
  set.seed(123L)
  edges <- rbind(
    c(1L, 2L), c(1L, 2L), c(1L, 3L), c(1L, 1L), c(2L, 3L),
    c(3L, 4L), c(3L, 4L), c(4L, 5L), c(4L, 6L), c(5L, 6L), c(5L, 5L),
    c(2L, 5L), c(1L, 6L)
  )
  b0 <- c(1L, 1L, 2L, 2L, 3L, 3L)
  w_list <- list(
    "real-exponential"   = c(1.2, 2.5, 0.7, 3.1, 1.5, 0.9, 2.1, 0.4, 1.8, 0.5, 2.2, 1.1, 0.8),
    "real-normal"        = c(-1.2, 2.5, 0.7, 3.1, -0.5, 0.9, 1.4, -2.1, 0.3, 1.1, -0.8, 0.2, -1.5),
    "discrete-geometric" = c(0, 2, 1, 4, 0, 3, 2, 1, 5, 1, 0, 2, 3),
    "discrete-poisson"   = c(1, 3, 0, 5, 2, 1, 4, 2, 0, 3, 1, 2, 0),
    "discrete-binomial"  = c(1, 3, 0, 4, 2, 1, 3, 2, 0, 4, 1, 2, 3)
  )

  for (dir in c(FALSE, TRUE)) {
    g <- graph_add_edges(graph_create(6, directed = dir), edges)
    for (tp in names(w_list)) {
      for (pm in c("microcanonical", "bayesian")) {
        st <- weighted_block_state(g, weights = w_list[[tp]], types = tp,
                                   params = pm, blocks = b0)
        errs <- rgraphtool:::cpp_weighted_check_moves(
          rgraphtool:::graph_pointer(g), st$blocks, st$weights, st$types,
          st$cpp_params, st$degree_corrected, TRUE, 40L, 1000
        )
        expect_lt(max(errs), 1e-9)
      }
    }
  }
})

test_that("weighted SBM recovers weight-only communities invisible to topology alone", {
  set.seed(42L)
  n <- 30L
  truth <- rep(1:2, each = n / 2L)
  # Uniform Erdős-Rényi topology (same edge probability p=0.45 within and between groups)
  ed_list <- list()
  w_norm <- numeric(0)
  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      if (stats::runif(1) < 0.45) {
        ed_list[[length(ed_list) + 1L]] <- c(i, j)
        mu <- if (truth[i] == truth[j]) 10.0 else -10.0
        w_norm <- c(w_norm, stats::rnorm(1, mean = mu, sd = 0.8))
      }
    }
  }
  edges <- do.call(rbind, ed_list)
  g <- graph_create(n, directed = FALSE)
  g <- graph_add_edges(g, edges, properties = list(w = w_norm))

  # Unweighted SBM should see no topological communities (B = 1)
  st_unw <- minimize_blockmodel_dl(g)
  expect_equal(length(block_sizes(st_unw)), 1L)

  # Flat and nested weighted SBMs should recover the 2 weight-defined communities cleanly
  st_w <- minimize_weighted_blockmodel_dl(g, weights = "w", types = "real-normal")
  expect_equal(st_w$B, 2L)
  expect_gt(partition_overlap(block_membership(st_w), truth, norm = TRUE), 0.95)
  expect_equal(sum(block_sizes(st_w)), n)
  expect_output(print(st_w), "weighted block state")

  nst_w <- minimize_nested_weighted_blockmodel_dl(g, weights = "w", types = "real-normal")
  expect_equal(nst_w$B[1L], 2L)
  expect_gt(partition_overlap(nested_levels(nst_w)[[1L]], truth, norm = TRUE), 0.95)
  expect_equal(length(nested_project(nst_w, level = 1L)), n)
  expect_s3_class(nested_summary(nst_w), "data.frame")
  expect_output(print(nst_w), "nested weighted block state")

  # B_min / B_max respected
  st_b3 <- minimize_weighted_blockmodel_dl(g, weights = "w", types = "real-normal",
                                           B_min = 3L, B_max = 3L)
  expect_equal(st_b3$B, 3L)
})

test_that("weighted MCMC, posterior sampling, and model_evidence work for flat and nested states", {
  set.seed(7L)
  n <- 20L
  truth <- rep(1:2, each = 10L)
  ed_list <- list()
  w_exp <- numeric(0)
  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      if (stats::runif(1) < 0.45) {
        ed_list[[length(ed_list) + 1L]] <- c(i, j)
        rate <- if (truth[i] == truth[j]) 0.1 else 2.0
        w_exp <- c(w_exp, stats::rexp(1, rate = rate))
      }
    }
  }
  g <- graph_add_edges(graph_create(n, directed = FALSE),
                       do.call(rbind, ed_list),
                       properties = list(w = w_exp))

  # Flat weighted MCMC + equilibration + annealing + evidence
  st <- minimize_weighted_blockmodel_dl(g, weights = "w", types = "real-exponential")
  sw <- mcmc_sweep(st, beta = 1.0, niter = 10L)
  expect_s3_class(sw$state, "rgraphtool_weighted_block_state")
  expect_equal(sw$entropy, weighted_description_length(sw$state), tolerance = 1e-9)

  eq <- mcmc_equilibrate(st, wait = 5L, nbreaks = 1L, max_niter = 15L)
  expect_s3_class(eq$state, "rgraphtool_weighted_block_state")

  an <- mcmc_anneal(st, beta_range = c(1, 5), niter = 3L,
                    mcmc_equilibrate_args = list(wait = 3L, max_niter = 6L))
  expect_s3_class(an$state, "rgraphtool_weighted_block_state")

  post <- sample_posterior(st, nsamples = 20L, burnin = 15L, thin = 2L)
  expect_equal(length(post$samples), 20L)
  expect_false(post$nested)

  ev <- model_evidence(st, nsamples = 20L, burnin = 15L, thin = 2L)
  expect_s3_class(ev, "rgraphtool_model_evidence")
  expect_true(is.finite(ev$log_evidence))
  expect_true(ev$posterior_entropy >= 0)

  # Nested weighted MCMC + posterior + evidence
  nst <- minimize_nested_weighted_blockmodel_dl(g, weights = "w", types = "real-exponential")
  nsw <- mcmc_sweep(nst, beta = 1.0, niter = 10L)
  expect_s3_class(nsw$state, "rgraphtool_nested_weighted_block_state")
  expect_equal(nsw$entropy, nested_weighted_description_length(nsw$state), tolerance = 1e-9)

  npost <- sample_posterior(nst, nsamples = 15L, burnin = 15L, thin = 2L)
  expect_equal(length(npost$samples), 15L)
  expect_true(npost$nested)

  nev <- model_evidence(nst, nsamples = 15L, burnin = 15L, thin = 2L)
  expect_s3_class(nev, "rgraphtool_model_evidence")
  expect_true(is.finite(nev$log_evidence))
})

test_that("compare_weight_models selects true weight model with Jacobian adjustment", {
  set.seed(99L)
  n <- 24L
  truth <- rep(1:2, each = 12L)
  ed_list <- list()
  w_ln <- numeric(0)
  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      if (stats::runif(1) < 0.5) {
        ed_list[[length(ed_list) + 1L]] <- c(i, j)
        mu <- if (truth[i] == truth[j]) 3.0 else 0.0
        w_ln <- c(w_ln, exp(stats::rnorm(1, mean = mu, sd = 0.9)))
      }
    }
  }
  g <- graph_add_edges(graph_create(n, directed = FALSE),
                       do.call(rbind, ed_list),
                       properties = list(w = w_ln))

  cmp <- compare_weight_models(g, weights = "w", nested = FALSE, method = "map")
  expect_s3_class(cmp, "rgraphtool_weight_model_comparison")
  expect_equal(cmp$best_model, "log-normal")
  expect_equal(sum(cmp$summary$posterior_weight), 1, tolerance = 1e-9)
  ln_row <- cmp$summary[cmp$summary$model == "log-normal", ]
  expect_equal(ln_row$jacobian_adj, sum(log(w_ln)), tolerance = 1e-10)
  expect_output(print(cmp), "Best weight model: log-normal")

  # Posterior evidence mode
  cmp_post <- compare_weight_models(
    g, weights = "w",
    candidates = c("real-exponential", "log-normal"),
    nested = FALSE, method = "posterior",
    nsamples = 20L, burnin = 15L, thin = 2L
  )
  expect_equal(cmp_post$best_model, "log-normal")
  expect_true(all(is.finite(cmp_post$summary$neg_log_evidence)))
})

test_that("weighted SBM validates domain, covariate specs, and edge cases", {
  g0 <- graph_create(4, directed = FALSE)
  st0 <- weighted_block_state(g0, weights = numeric(0), types = "real-normal")
  expect_true(is.finite(weighted_description_length(st0)))

  g <- graph_add_edges(graph_create(4, directed = FALSE),
                       rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L)))
  expect_error(weighted_block_state(g, weights = c(-1, 2, 3), types = "real-exponential"),
               "non-negative")
  expect_error(weighted_block_state(g, weights = c(1.5, 2, 3), types = "discrete-poisson"),
               "integer")
  expect_error(weighted_block_state(g, weights = c(1, 2, 5), types = "discrete-binomial",
                                    params = list(N = 3)),
               "N")
  expect_error(weighted_block_state(g, weights = c(1, 2), types = "real-normal"),
               "length")
  expect_error(weighted_block_state(g, weights = c(1, NA, 3), types = "real-normal"),
               "finite|missing")
  expect_error(weighted_block_state(g, weights = c(1, 2, 3), types = "unknown-dist"),
               "unknown covariate type")

  # Mutating the graph invalidates previously created weighted states
  st <- weighted_block_state(g, weights = c(1, 2, 3), types = "real-normal")
  graph_add_edges(g, rbind(c(1L, 4L)))
  expect_error(weighted_description_length(st), "changed")
})
