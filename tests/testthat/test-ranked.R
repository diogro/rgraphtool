test_that("flat and nested RankedBlockState matches pinned graph-tool 3.9 values", {
  edges <- rbind(
    c(1L, 2L), c(1L, 2L), c(2L, 1L), c(1L, 1L),
    c(1L, 3L), c(2L, 3L), c(2L, 4L), c(3L, 2L),
    c(3L, 4L), c(3L, 5L), c(4L, 5L), c(4L, 6L), c(5L, 6L), c(6L, 5L), c(6L, 4L),
    c(1L, 5L), c(5L, 1L), c(6L, 6L)
  )
  b <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L)
  bs <- list(b, c(1L, 1L, 2L, 2L), c(1L, 1L))
  g <- graph_add_edges(graph_create(7L, directed = TRUE), edges)

  ref <- list(
    list(u = c(0, 0, 2, 1), dist = FALSE, uniq = FALSE, dc = TRUE,
         flat = 56.823150531968, nolbl = 51.476043001251, nest = 62.420457041883,
         lv = c(41.462884593076, 15.017401196086, 5.940171252720)),
    list(u = c(0, 0, 2, 1), dist = FALSE, uniq = FALSE, dc = FALSE,
         flat = 53.148662323888, nolbl = 47.801554793170, nest = 58.745968833802,
         lv = c(37.788396384996, 15.017401196086, 5.940171252720)),
    list(u = c(0, 0, 2, 1), dist = TRUE,  uniq = FALSE, dc = TRUE,
         flat = 60.288886434768, nolbl = 54.941778904050, nest = 65.886192944682,
         lv = c(44.928620495876, 15.017401196086, 5.940171252720)),
    list(u = c(0, 0, 2, 1), dist = TRUE,  uniq = FALSE, dc = FALSE,
         flat = 56.614398226687, nolbl = 51.267290695970, nest = 62.211704736602,
         lv = c(41.254132287796, 15.017401196086, 5.940171252720)),
    list(u = c(0, 1, 3, 2), dist = FALSE, uniq = FALSE, dc = TRUE,
         flat = 56.581988475151, nolbl = 51.234880944434, nest = 62.179294985066,
         lv = c(41.221722536259, 15.017401196086, 5.940171252720)),
    list(u = c(0, 1, 3, 2), dist = FALSE, uniq = FALSE, dc = FALSE,
         flat = 52.907500267071, nolbl = 47.560392736353, nest = 58.504806776985,
         lv = c(37.547234328179, 15.017401196086, 5.940171252720)),
    list(u = c(0, 1, 3, 2), dist = FALSE, uniq = TRUE,  dc = TRUE,
         flat = 51.234880944434, nolbl = 51.234880944434, nest = 56.832187454348,
         lv = c(35.874615005542, 15.017401196086, 5.940171252720)),
    list(u = c(0, 1, 3, 2), dist = FALSE, uniq = TRUE,  dc = FALSE,
         flat = 47.560392736353, nolbl = 47.560392736353, nest = 53.157699246268,
         lv = c(32.200126797461, 15.017401196086, 5.940171252720)),
    list(u = c(0, 1, 3, 2), dist = TRUE,  uniq = FALSE, dc = TRUE,
         flat = 60.406669470424, nolbl = 55.059561939707, nest = 66.003975980339,
         lv = c(45.046403531532, 15.017401196086, 5.940171252720)),
    list(u = c(0, 1, 3, 2), dist = TRUE,  uniq = FALSE, dc = FALSE,
         flat = 56.732181262344, nolbl = 51.385073731626, nest = 62.329487772258,
         lv = c(41.371915323452, 15.017401196086, 5.940171252720)),
    list(u = c(0, 1, 3, 2), dist = TRUE,  uniq = TRUE,  dc = TRUE,
         flat = 54.877240382913, nolbl = 54.877240382913, nest = 60.474546892827,
         lv = c(39.516974444021, 15.017401196086, 5.940171252720)),
    list(u = c(0, 1, 3, 2), dist = TRUE,  uniq = TRUE,  dc = FALSE,
         flat = 51.202752174832, nolbl = 51.202752174832, nest = 56.800058684747,
         lv = c(35.842486235940, 15.017401196086, 5.940171252720))
  )

  for (r in ref) {
    st <- ranked_block_state(g, blocks = b, u = r$u, distance = r$dist,
                             unique = r$uniq, degree_corrected = r$dc)
    expect_equal(ranked_description_length(st, label_dl = TRUE), r$flat, tolerance = 1e-10)
    expect_equal(ranked_description_length(st, label_dl = FALSE), r$nolbl, tolerance = 1e-10)

    nst <- nested_ranked_block_state(g, levels = bs, u = r$u, distance = r$dist,
                                     unique = r$uniq, degree_corrected = r$dc)
    expect_equal(nested_ranked_description_length(nst), r$nest, tolerance = 1e-10)
    expect_equal(nested_ranked_level_entropy(nst), r$lv, tolerance = 1e-10)

    # Incremental vertex moves AND rank relabels match full recomputation
    set.seed(11L)
    errs <- graphtool:::cpp_ranked_check_moves(
      graphtool:::graph_pointer(g), st$blocks, st$u, st$distance, st$unique,
      as.numeric(st$M), st$degree_corrected, TRUE, 36L, 1000
    )
    expect_lt(max(errs), 1e-9)
  }
})

test_that("ranked SBM recovers ordered feed-forward hierarchy and beats unordered SBM", {
  set.seed(42L)
  # 6-tier ordered communities: cohesive within-tier modules + strong upstream flow
  # from lower tiers to higher tiers (Peixoto 2022 ordered community model)
  n <- 72L
  tier <- rep(1:6, each = 12L)
  ed <- list()
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      ti <- tier[i]; tj <- tier[j]
      p <- if (ti == tj) 0.70 else if (tj > ti) 0.20 else 0.005
      if (stats::runif(1) < p) ed[[length(ed) + 1L]] <- c(i, j)
    }
  }
  g <- graph_add_edges(graph_create(n, directed = TRUE), do.call(rbind, ed))

  st <- minimize_ranked_blockmodel_dl(g)
  expect_gt(partition_overlap(block_membership(st), tier, norm = TRUE), 0.90)
  expect_equal(sum(block_sizes(st)), n)
  expect_output(print(st), "ranked block state")

  # Inferred vertex rank correlates strongly with true tier ordering
  ord_v <- ranked_vertex_order(st, normalize = TRUE)
  expect_gt(abs(stats::cor(ord_v, tier)), 0.90)
  ec <- ranked_edge_counts(st)
  expect_gt(ec$upstream, 5 * ec$downstream)
  dirs <- ranked_edge_directions(st)
  expect_equal(sum(dirs == 1L), ec$upstream)
  expect_equal(sum(dirs == 0L), ec$lateral)
  expect_equal(sum(dirs == -1L), ec$downstream)

  nst <- minimize_nested_ranked_blockmodel_dl(g)
  expect_gt(partition_overlap(nested_levels(nst)[[1L]], tier, norm = TRUE), 0.90)
  expect_equal(length(nested_project(nst, 0L)), n)
  expect_s3_class(nested_summary(nst), "data.frame")
  expect_output(print(nst), "nested ranked block state")

  # Ordered model beats unordered SBM on an ordered community hierarchy
  cmp <- compare_ordered_models(g, method = "map")
  expect_s3_class(cmp, "graphtool_ordered_model_comparison")
  expect_true(cmp$best_model %in% c("ranked_flat", "ranked_nested"))
  expect_equal(sum(cmp$summary$posterior_weight), 1, tolerance = 1e-9)
  expect_output(print(cmp), "Best model:")
})

test_that("ranked MCMC, posterior sampling, model_evidence, and ranked layouts work", {
  set.seed(19L)
  n <- 18L
  tier <- rep(1:3, each = 6L)
  ed <- list()
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      p <- if (tier[j] > tier[i]) 0.65 else if (tier[i] == tier[j]) 0.15 else 0.02
      if (stats::runif(1) < p) ed[[length(ed) + 1L]] <- c(i, j)
    }
  }
  g <- graph_add_edges(graph_create(n, directed = TRUE), do.call(rbind, ed))

  st <- minimize_ranked_blockmodel_dl(g)
  sw <- mcmc_sweep(st, beta = 1.0, niter = 10L)
  expect_s3_class(sw$state, "graphtool_ranked_block_state")
  expect_equal(sw$entropy, ranked_description_length(sw$state), tolerance = 1e-9)

  eq <- mcmc_equilibrate(st, wait = 5L, nbreaks = 1L, max_niter = 15L)
  expect_s3_class(eq$state, "graphtool_ranked_block_state")

  an <- mcmc_anneal(st, beta_range = c(1, 5), niter = 3L,
                    mcmc_equilibrate_args = list(wait = 3L, max_niter = 6L))
  expect_s3_class(an$state, "graphtool_ranked_block_state")

  post <- sample_posterior(st, nsamples = 20L, burnin = 15L, thin = 2L)
  expect_equal(length(post$samples), 20L)
  expect_false(post$nested)

  ev <- model_evidence(st, nsamples = 20L, burnin = 15L, thin = 2L)
  expect_s3_class(ev, "graphtool_model_evidence")
  expect_true(is.finite(ev$log_evidence))

  nst <- minimize_nested_ranked_blockmodel_dl(g)
  nsw <- mcmc_sweep(nst, beta = 1.0, niter = 10L)
  expect_s3_class(nsw$state, "graphtool_nested_ranked_block_state")
  expect_equal(nsw$entropy, nested_ranked_description_length(nsw$state), tolerance = 1e-9)

  nev <- model_evidence(nst, nsamples = 15L, burnin = 15L, thin = 2L)
  expect_s3_class(nev, "graphtool_model_evidence")
  expect_true(is.finite(nev$log_evidence))

  cmp_p <- compare_ordered_models(g, models = c("unordered_flat", "ranked_flat"),
                                  method = "posterior", nsamples = 20L,
                                  burnin = 15L, thin = 2L)
  expect_true(all(is.finite(cmp_p$summary$neg_log_evidence)))

  # Both plot_ranked_hierarchy (flat & nested) and plot_hierarchy (nested) succeed
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
  expect_silent(plot_ranked_hierarchy(st, main = "Flat ranked layout"))
  expect_silent(plot_ranked_hierarchy(nst, main = "Nested ranked layout"))
  expect_silent(plot_hierarchy(nst, main = "Nested circular hierarchy"))
})

test_that("ranked SBM validates directedness, rank constraints, and edge cases", {
  g_u <- graph_add_edges(graph_create(4, directed = FALSE), rbind(c(1L, 2L), c(2L, 3L)))
  expect_error(ranked_block_state(g_u), "directed")
  expect_error(minimize_ranked_blockmodel_dl(g_u), "directed")
  expect_error(nested_ranked_block_state(g_u), "directed")
  expect_error(compare_ordered_models(g_u), "directed")

  g_d <- graph_add_edges(graph_create(4, directed = TRUE), rbind(c(1L, 2L), c(2L, 3L)))
  expect_error(ranked_block_state(g_d, blocks = c(1L, 1L, 2L, 2L), u = c(0, 0), unique = TRUE),
               "distinct")
  expect_error(ranked_block_state(g_d, blocks = c(1L, 1L, 2L, 2L), u = c(0)),
               "length")
  expect_error(compare_ordered_models(g_d, models = "bogus"), "models must be")

  st <- ranked_block_state(g_d, blocks = c(1L, 1L, 2L, 2L))
  graph_add_edges(g_d, rbind(c(3L, 4L)))
  expect_error(ranked_description_length(st), "changed")
})
