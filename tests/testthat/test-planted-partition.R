test_that("PPBlockState matches pinned graph-tool 3.9 values across all 16 configurations", {
  edges <- rbind(
    c(1L, 2L), c(1L, 2L), c(1L, 3L), c(1L, 1L), c(2L, 3L),
    c(3L, 4L), c(4L, 5L), c(4L, 5L), c(5L, 6L), c(6L, 4L), c(5L, 5L)
  )
  b <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L)

  ref <- list(
    list(dir = FALSE, unif = FALSE, dc = TRUE,  dk = "distributed", S = 30.045371343521),
    list(dir = FALSE, unif = FALSE, dc = TRUE,  dk = "uniform",     S = 30.323239794524),
    list(dir = FALSE, unif = FALSE, dc = FALSE, dk = "distributed", S = 27.804380235024),
    list(dir = FALSE, unif = FALSE, dc = FALSE, dk = "uniform",     S = 27.804380235024),
    list(dir = FALSE, unif = TRUE,  dc = TRUE,  dk = "distributed", S = 31.312410400664),
    list(dir = FALSE, unif = TRUE,  dc = TRUE,  dk = "uniform",     S = 31.590278851668),
    list(dir = FALSE, unif = TRUE,  dc = FALSE, dk = "distributed", S = 29.071419292167),
    list(dir = FALSE, unif = TRUE,  dc = FALSE, dk = "uniform",     S = 29.071419292167),
    list(dir = TRUE,  unif = FALSE, dc = TRUE,  dk = "distributed", S = 36.369463674785),
    list(dir = TRUE,  unif = FALSE, dc = TRUE,  dk = "uniform",     S = 36.924625112371),
    list(dir = TRUE,  unif = FALSE, dc = FALSE, dk = "distributed", S = 34.042704860064),
    list(dir = TRUE,  unif = FALSE, dc = FALSE, dk = "uniform",     S = 34.042704860064),
    list(dir = TRUE,  unif = TRUE,  dc = TRUE,  dk = "distributed", S = 37.636502731928),
    list(dir = TRUE,  unif = TRUE,  dc = TRUE,  dk = "uniform",     S = 38.191664169514),
    list(dir = TRUE,  unif = TRUE,  dc = FALSE, dk = "distributed", S = 35.309743917207),
    list(dir = TRUE,  unif = TRUE,  dc = FALSE, dk = "uniform",     S = 35.309743917207)
  )

  g_u <- graph_add_edges(graph_create(7, directed = FALSE), edges)
  g_d <- graph_add_edges(graph_create(7, directed = TRUE), edges)

  for (r in ref) {
    g <- if (r$dir) g_d else g_u
    st <- pp_block_state(g, blocks = b, uniform = r$unif,
                         degree_corrected = r$dc, degree_prior = r$dk)
    expect_equal(pp_description_length(st), r$S, tolerance = 1e-10)
    set.seed(11L)
    errs <- graphtool:::cpp_pp_check_moves(
      graphtool:::graph_pointer(g), st$blocks, st$uniform, st$degree_corrected,
      r$dk == "distributed", 30L, 1000
    )
    expect_lt(max(errs), 1e-9)
  }
})

test_that("planted-partition recovers assortative communities and rejects disassortative ones", {
  set.seed(42L)
  # 1. Strongly assortative 3-community graph
  n <- 36L
  truth <- rep(1:3, each = 12L)
  ed_a <- list()
  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      p <- if (truth[i] == truth[j]) 0.70 else 0.03
      if (stats::runif(1) < p) ed_a[[length(ed_a) + 1L]] <- c(i, j)
    }
  }
  g_assort <- graph_add_edges(graph_create(n, directed = FALSE), do.call(rbind, ed_a))

  st_pp <- minimize_pp_blockmodel_dl(g_assort)
  expect_equal(st_pp$B, 3L)
  expect_gt(partition_overlap(block_membership(st_pp), truth, norm = TRUE), 0.95)
  expect_equal(sum(block_sizes(st_pp)), n)
  asrt <- pp_assortativity(st_pp)
  expect_gt(asrt$fraction_in, 0.85)
  expect_equal(asrt$E_in + asrt$E_out, asrt$E_total)
  expect_output(print(st_pp), "planted-partition block state")

  # Assortative model comparison favors PP on a pure planted-partition graph
  cmp_a <- compare_assortative_models(g_assort, method = "map")
  expect_s3_class(cmp_a, "graphtool_assortative_model_comparison")
  expect_true(cmp_a$best_model %in% c("pp", "pp_uniform"))
  expect_equal(sum(cmp_a$summary$posterior_weight), 1, tolerance = 1e-9)
  expect_output(print(cmp_a), "Best model:")

  # 2. Mixed assortative + bipartite 4-block network:
  #    blocks 1 & 2 are dense assortative cliques, while blocks 3 & 4 form a
  #    dense bipartite pair. No single (E_in, E_out) rate can explain both, so
  #    the general SBM decisively beats PP.
  ed_b <- list()
  truth_b <- rep(1:4, each = 10L)
  P_mat <- rbind(
    c(0.75, 0.02, 0.02, 0.02),
    c(0.02, 0.75, 0.02, 0.02),
    c(0.02, 0.02, 0.02, 0.75),
    c(0.02, 0.02, 0.75, 0.02)
  )
  for (i in 1:39) {
    for (j in (i + 1L):40) {
      if (stats::runif(1) < P_mat[truth_b[i], truth_b[j]]) {
        ed_b[[length(ed_b) + 1L]] <- c(i, j)
      }
    }
  }
  g_nonassort <- graph_add_edges(graph_create(40L, directed = FALSE), do.call(rbind, ed_b))
  cmp_b <- compare_assortative_models(g_nonassort, method = "map")
  expect_true(cmp_b$best_model %in% c("sbm_flat", "sbm_nested"))
})

test_that("PPBlockState supports MCMC, posterior sampling, model_evidence, and validation", {
  set.seed(9L)
  edges <- rbind(t(utils::combn(1:8, 2)), t(utils::combn(9:16, 2)), c(8, 9), c(1, 10))
  g <- graph_add_edges(graph_create(16L, directed = FALSE), edges)
  st <- minimize_pp_blockmodel_dl(g)

  sw <- mcmc_sweep(st, beta = 1.0, niter = 10L)
  expect_s3_class(sw$state, "graphtool_pp_block_state")
  expect_equal(sw$entropy, pp_description_length(sw$state), tolerance = 1e-9)

  eq <- mcmc_equilibrate(st, wait = 5L, nbreaks = 1L, max_niter = 15L)
  expect_s3_class(eq$state, "graphtool_pp_block_state")

  an <- mcmc_anneal(st, beta_range = c(1, 5), niter = 3L,
                    mcmc_equilibrate_args = list(wait = 3L, max_niter = 6L))
  expect_s3_class(an$state, "graphtool_pp_block_state")

  post <- sample_posterior(st, nsamples = 20L, burnin = 15L, thin = 2L)
  expect_equal(length(post$samples), 20L)
  expect_false(post$nested)

  ev <- model_evidence(st, nsamples = 20L, burnin = 15L, thin = 2L)
  expect_s3_class(ev, "graphtool_model_evidence")
  expect_true(is.finite(ev$log_evidence))

  cmp_p <- compare_assortative_models(g, models = c("pp", "sbm_flat"),
                                      method = "posterior", nsamples = 20L,
                                      burnin = 15L, thin = 2L)
  expect_true(all(is.finite(cmp_p$summary$neg_log_evidence)))

  # Validation & graph mutation
  g0 <- graph_create(4, directed = FALSE)
  expect_true(is.finite(pp_description_length(pp_block_state(g0))))
  expect_error(pp_block_state(g, uniform = NA), "TRUE or FALSE")
  expect_error(compare_assortative_models(g, models = "invalid"), "models must be")
  graph_add_edges(g, rbind(c(2L, 12L)))
  expect_error(pp_description_length(st), "changed")
})
