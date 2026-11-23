test_that("flat and nested layered SBM matches pinned graph-tool 3.9 values across all 16 configurations", {
  # 7 vertices (vertex 7 isolated), 3 layers, parallel edges within & across layers,
  # self-loops, and block 3 absent from layer 2
  edges <- rbind(
    c(1L, 2L), c(1L, 2L), c(1L, 2L), c(1L, 1L),
    c(2L, 3L), c(3L, 4L), c(3L, 4L), c(4L, 5L),
    c(4L, 6L), c(5L, 6L), c(5L, 5L), c(2L, 5L)
  )
  ec <- c(1L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 3L, 3L, 3L, 1L)
  b_flat <- c(1L, 1L, 2L, 2L, 3L, 3L, 3L)
  bs_nest <- list(b_flat, c(1L, 1L, 2L), c(1L, 1L))

  ref <- list(
    list(dir = FALSE, indep = FALSE, dc = TRUE,  dk = "distributed",
         flat = 42.115728394612, nolm = 42.115728394612, nest = 46.252748520003,
         lv = c(23.922126936827, 9.692766521220, 12.637855061956)),
    list(dir = FALSE, indep = FALSE, dc = TRUE,  dk = "uniform",
         flat = 42.434182125730, nolm = 42.434182125730, nest = 46.571202251122,
         lv = c(24.240580667946, 9.692766521220, 12.637855061956)),
    list(dir = FALSE, indep = FALSE, dc = FALSE, dk = "distributed",
         flat = 41.717814876268, nolm = 41.717814876268, nest = 45.854835001659,
         lv = c(23.524213418483, 9.692766521220, 12.637855061956)),
    list(dir = FALSE, indep = FALSE, dc = FALSE, dk = "uniform",
         flat = 41.717814876268, nolm = 41.717814876268, nest = 45.854835001659,
         lv = c(23.524213418483, 9.692766521220, 12.637855061956)),
    list(dir = FALSE, indep = TRUE,  dc = TRUE,  dk = "distributed",
         flat = 40.161850075907, nolm = 37.859264982913, nest = 41.303137927745,
         lv = c(21.968248618123, 8.999619340661, 10.335269968962)),
    list(dir = FALSE, indep = TRUE,  dc = TRUE,  dk = "uniform",
         flat = 42.071392580792, nolm = 39.768807487798, nest = 43.212680432629,
         lv = c(23.877791123007, 8.999619340661, 10.335269968962)),
    list(dir = FALSE, indep = TRUE,  dc = FALSE, dk = "distributed",
         flat = 39.102702670825, nolm = 36.800117577831, nest = 40.243990522662,
         lv = c(20.909101213040, 8.999619340661, 10.335269968962)),
    list(dir = FALSE, indep = TRUE,  dc = FALSE, dk = "uniform",
         flat = 39.102702670825, nolm = 36.800117577831, nest = 40.243990522662,
         lv = c(20.909101213040, 8.999619340661, 10.335269968962)),
    list(dir = TRUE,  indep = FALSE, dc = TRUE,  dk = "distributed",
         flat = 51.272265537362, nolm = 51.272265537362, nest = 55.401647610257,
         lv = c(29.045674336345, 11.233211562168, 15.122761711744)),
    list(dir = TRUE,  indep = FALSE, dc = TRUE,  dk = "uniform",
         flat = 51.110563267771, nolm = 51.110563267771, nest = 55.239945340665,
         lv = c(28.883972066754, 11.233211562168, 15.122761711744)),
    list(dir = TRUE,  indep = FALSE, dc = FALSE, dk = "distributed",
         flat = 49.909687702860, nolm = 49.909687702860, nest = 54.039069775754,
         lv = c(27.683096501843, 11.233211562168, 15.122761711744)),
    list(dir = TRUE,  indep = FALSE, dc = FALSE, dk = "uniform",
         flat = 49.909687702860, nolm = 49.909687702860, nest = 54.039069775754,
         lv = c(27.683096501843, 11.233211562168, 15.122761711744)),
    list(dir = TRUE,  indep = TRUE,  dc = TRUE,  dk = "distributed",
         flat = 48.734114373055, nolm = 46.431529280061, nest = 49.174616991836,
         lv = c(26.507523172039, 10.540064381608, 12.127029438190)),
    list(dir = TRUE,  indep = TRUE,  dc = TRUE,  dk = "uniform",
         flat = 47.871068155700, nolm = 45.568483062706, nest = 48.311570774481,
         lv = c(25.644476954683, 10.540064381608, 12.127029438190)),
    list(dir = TRUE,  indep = TRUE,  dc = FALSE, dk = "distributed",
         flat = 47.294575497417, nolm = 44.991990404423, nest = 47.735078116197,
         lv = c(25.067984296400, 10.540064381608, 12.127029438190)),
    list(dir = TRUE,  indep = TRUE,  dc = FALSE, dk = "uniform",
         flat = 47.294575497417, nolm = 44.991990404423, nest = 47.735078116197,
         lv = c(25.067984296400, 10.540064381608, 12.127029438190))
  )

  g_u <- graph_add_edges(graph_create(7, directed = FALSE), edges,
                         properties = list(layer = ec))
  g_d <- graph_add_edges(graph_create(7, directed = TRUE), edges,
                         properties = list(layer = ec))

  for (r in ref) {
    g <- if (r$dir) g_d else g_u
    st <- layered_block_state(g, layers = "layer", independent = r$indep,
                              blocks = b_flat, degree_corrected = r$dc,
                              degree_prior = r$dk)
    expect_equal(layered_description_length(st, lmemb = TRUE), r$flat,
                 tolerance = 1e-10)
    expect_equal(layered_description_length(st, lmemb = FALSE), r$nolm,
                 tolerance = 1e-10)

    nst <- nested_layered_block_state(g, layers = "layer", independent = r$indep,
                                      levels = bs_nest, degree_corrected = r$dc,
                                      degree_prior = r$dk)
    expect_equal(nested_layered_description_length(nst), r$nest, tolerance = 1e-10)
    expect_equal(nested_layered_level_entropy(nst), r$lv, tolerance = 1e-10)

    # Incremental move deltas match full recomputation
    set.seed(11L)
    errs <- graphtool:::cpp_layered_check_moves(
      graphtool:::graph_pointer(g), st$blocks, st$layers, st$independent,
      st$degree_corrected, r$dk == "distributed", 30L, 1000
    )
    expect_lt(max(errs), 1e-9)
  }
})

test_that("layered_block_state accepts character/factor layers and netzschleuder string layers", {
  edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 1L), c(1L, 3L))
  g <- graph_add_edges(graph_create(3, directed = FALSE), edges)
  lbls <- c("marriage", "business", "marriage", "business")
  st_chr <- layered_block_state(g, layers = lbls, independent = TRUE,
                                blocks = c(1L, 1L, 2L))
  st_fac <- layered_block_state(g, layers = factor(lbls), independent = TRUE,
                                blocks = c(1L, 1L, 2L))
  expect_equal(st_chr$C, 2L)
  expect_setequal(st_chr$layer_levels, c("business", "marriage"))
  expect_equal(layered_description_length(st_chr),
               layered_description_length(st_fac), tolerance = 1e-12)

  # Single-layer graph matches unweighted SBM description length exactly
  st_one <- layered_block_state(g, layers = rep("only", 4L), independent = FALSE,
                                blocks = c(1L, 1L, 2L))
  st_unw <- block_state(g, blocks = c(1L, 1L, 2L), degree_corrected = TRUE)
  expect_equal(layered_description_length(st_one),
               block_description_length(st_unw), tolerance = 1e-10)

  # Netzschleuder archive with categorical string layer column
  td <- tempfile("ns_layered_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  writeLines(c("# source, target, layer",
               "0, 1, marriage",
               "1, 2, business",
               "2, 0, marriage"), file.path(td, "edges.csv"))
  writeLines(c("# index, name", "0, A", "1, B", "2, C"), file.path(td, "nodes.csv"))
  zfile <- file.path(td, "toy_layered.csv.zip")
  oldwd <- setwd(td)
  on.exit(setwd(oldwd), add = TRUE)
  utils::zip(zfile, files = c("edges.csv", "nodes.csv"), flags = "-q")
  setwd(oldwd)

  net <- netzschleuder_network("toy_layered", file = zfile)
  expect_true("layer" %in% graph_edge_properties(net$graph))
  expect_setequal(net$edge_levels$layer, c("business", "marriage"))
  st_net <- layered_block_state(net$graph, layers = "layer", independent = TRUE)
  expect_equal(st_net$C, 2L)
})

test_that("layered SBM recovers communities invisible in the collapsed network", {
  set.seed(42L)
  n <- 30L
  truth <- rep(1:2, each = n / 2L)
  # Layer 1 ("friend"): edges predominantly within groups (p_in = 0.55, p_out = 0.05)
  # Layer 2 ("foe"):    edges predominantly between groups (p_in = 0.05, p_out = 0.55)
  # Collapsed graph has uniform edge probability 0.60 everywhere -> invisible to unweighted SBM!
  ed_list <- list()
  lay <- character(0)
  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      same <- (truth[i] == truth[j])
      if (stats::runif(1) < (if (same) 0.55 else 0.05)) {
        ed_list[[length(ed_list) + 1L]] <- c(i, j)
        lay <- c(lay, "friend")
      }
      if (stats::runif(1) < (if (same) 0.05 else 0.55)) {
        ed_list[[length(ed_list) + 1L]] <- c(i, j)
        lay <- c(lay, "foe")
      }
    }
  }
  g <- graph_add_edges(graph_create(n, directed = FALSE), do.call(rbind, ed_list))

  st_unw <- minimize_blockmodel_dl(g)
  expect_equal(length(block_sizes(st_unw)), 1L)

  for (indep in c(FALSE, TRUE)) {
    st_l <- minimize_layered_blockmodel_dl(g, layers = lay, independent = indep)
    expect_equal(st_l$B, 2L)
    expect_gt(partition_overlap(block_membership(st_l), truth, norm = TRUE), 0.95)
    expect_equal(sum(block_sizes(st_l)), n)
    expect_output(print(st_l), "layered block state")

    nst_l <- minimize_nested_layered_blockmodel_dl(g, layers = lay, independent = indep)
    expect_equal(nst_l$B[1L], 2L)
    expect_gt(partition_overlap(nested_levels(nst_l)[[1L]], truth, norm = TRUE), 0.95)
    expect_equal(length(nested_project(nst_l, level = 0L)), n)
    expect_s3_class(nested_summary(nst_l), "data.frame")
    expect_output(print(nst_l), "nested layered block state")
  }

  # B_min / B_max bounds
  st_b3 <- minimize_layered_blockmodel_dl(g, layers = lay, independent = FALSE,
                                          B_min = 3L, B_max = 3L)
  expect_equal(st_b3$B, 3L)
})

test_that("layered MCMC, posterior sampling, model_evidence, and plot_hierarchy work", {
  set.seed(17L)
  n <- 18L
  truth <- rep(1:2, each = 9L)
  ed_list <- list()
  lay <- integer(0)
  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      same <- (truth[i] == truth[j])
      if (stats::runif(1) < (if (same) 0.55 else 0.05)) {
        ed_list[[length(ed_list) + 1L]] <- c(i, j)
        lay <- c(lay, 1L)
      }
      if (stats::runif(1) < (if (same) 0.05 else 0.55)) {
        ed_list[[length(ed_list) + 1L]] <- c(i, j)
        lay <- c(lay, 2L)
      }
    }
  }
  g <- graph_add_edges(graph_create(n, directed = FALSE), do.call(rbind, ed_list),
                       properties = list(layer = lay))

  st <- minimize_layered_blockmodel_dl(g, layers = "layer", independent = TRUE)
  sw <- mcmc_sweep(st, beta = 1.0, niter = 10L)
  expect_s3_class(sw$state, "graphtool_layered_block_state")
  expect_equal(sw$entropy, layered_description_length(sw$state), tolerance = 1e-9)

  eq <- mcmc_equilibrate(st, wait = 5L, nbreaks = 1L, max_niter = 15L)
  expect_s3_class(eq$state, "graphtool_layered_block_state")

  an <- mcmc_anneal(st, beta_range = c(1, 5), niter = 3L,
                    mcmc_equilibrate_args = list(wait = 3L, max_niter = 6L))
  expect_s3_class(an$state, "graphtool_layered_block_state")

  post <- sample_posterior(st, nsamples = 20L, burnin = 15L, thin = 2L)
  expect_equal(length(post$samples), 20L)
  expect_false(post$nested)

  ev <- model_evidence(st, nsamples = 20L, burnin = 15L, thin = 2L)
  expect_s3_class(ev, "graphtool_model_evidence")
  expect_true(is.finite(ev$log_evidence))
  expect_true(ev$posterior_entropy >= 0)

  # Nested layered MCMC + posterior + evidence + plot_hierarchy
  nst <- minimize_nested_layered_blockmodel_dl(g, layers = "layer", independent = TRUE)
  nsw <- mcmc_sweep(nst, beta = 1.0, niter = 10L)
  expect_s3_class(nsw$state, "graphtool_nested_layered_block_state")
  expect_equal(nsw$entropy, nested_layered_description_length(nsw$state), tolerance = 1e-9)

  npost <- sample_posterior(nst, nsamples = 15L, burnin = 15L, thin = 2L)
  expect_equal(length(npost$samples), 15L)
  expect_true(npost$nested)

  nev <- model_evidence(nst, nsamples = 15L, burnin = 15L, thin = 2L)
  expect_s3_class(nev, "graphtool_model_evidence")
  expect_true(is.finite(nev$log_evidence))

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
  expect_silent(plot_hierarchy(nst, main = "Test layered hierarchy"))
})

test_that("compare_layered_models distinguishes independent layers from edge covariates", {
  set.seed(101L)
  # Generative independent-layers network: 4 sparse layers on disjoint subsets of vertices
  n <- 32L
  ed_list <- list()
  lay <- integer(0)
  for (l in 1:4) {
    verts <- ((l - 1L) * 8L + 1L):(l * 8L)
    for (i in seq_len(length(verts) - 1L)) {
      for (j in (i + 1L):length(verts)) {
        if (stats::runif(1) < 0.75) {
          ed_list[[length(ed_list) + 1L]] <- c(verts[i], verts[j])
          lay <- c(lay, l)
        }
      }
    }
  }
  g_indep <- graph_add_edges(graph_create(n, directed = FALSE), do.call(rbind, ed_list),
                             properties = list(layer = lay))

  cmp_map <- compare_layered_models(g_indep, layers = "layer", nested = TRUE,
                                    method = "map")
  expect_s3_class(cmp_map, "graphtool_layered_model_comparison")
  expect_equal(cmp_map$best_model, "independent_layers")
  expect_equal(sum(cmp_map$summary$posterior_weight), 1, tolerance = 1e-9)
  expect_output(print(cmp_map), "Best layered model: independent_layers")

  cmp_post <- compare_layered_models(g_indep, layers = "layer", nested = FALSE,
                                     method = "posterior", nsamples = 20L,
                                     burnin = 15L, thin = 2L)
  expect_equal(cmp_post$best_model, "independent_layers")
  expect_true(all(is.finite(cmp_post$summary$neg_log_evidence)))

  # Generative edge-covariates network: shared degree sequence with all vertices active
  # across 2 layers (random categorical edge labels drawn independently of degrees)
  n2 <- 30L
  ed2 <- list()
  lay2 <- integer(0)
  for (i in 1:(n2 - 1L)) {
    for (j in (i + 1L):n2) {
      if (stats::runif(1) < 0.75) {
        ed2[[length(ed2) + 1L]] <- c(i, j)
        lay2 <- c(lay2, sample.int(2L, 1L))
      }
    }
  }
  g_cov <- graph_add_edges(graph_create(n2, directed = FALSE), do.call(rbind, ed2),
                           properties = list(layer = lay2))
  cmp_cov <- compare_layered_models(g_cov, layers = "layer", nested = FALSE,
                                    method = "map")
  expect_equal(cmp_cov$best_model, "edge_covariates")
})

test_that("layered SBM validates inputs and edge cases", {
  g0 <- graph_create(4, directed = FALSE)
  st0 <- layered_block_state(g0, layers = integer(0), independent = TRUE)
  expect_true(is.finite(layered_description_length(st0)))

  g <- graph_add_edges(graph_create(4, directed = FALSE),
                       rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L)))
  expect_error(layered_block_state(g, layers = c(1L, 2L)), "length")
  expect_error(layered_block_state(g, layers = c(1L, NA, 2L)), "missing")
  expect_error(layered_block_state(g, layers = c(1L, 2L, 1L), independent = NA),
               "TRUE or FALSE")
  expect_error(layered_block_state(g, layers = "nonexistent"), "does not exist")
  expect_error(compare_layered_models(g, layers = c(1L, 2L, 1L), models = "bogus"),
               "models must be")

  # Graph mutation invalidates existing layered state
  st <- layered_block_state(g, layers = c(1L, 2L, 1L))
  graph_add_edges(g, rbind(c(1L, 4L)))
  expect_error(layered_description_length(st), "changed")
})
