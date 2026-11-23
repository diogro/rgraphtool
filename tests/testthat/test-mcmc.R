test_that("exact MCMC proposal probabilities normalize and match post-move reverse probs", {
  set.seed(101)
  for (directed in c(FALSE, TRUE)) {
    g <- graph_create(20, directed = directed)
    graph_add_edges(g, rbind(cbind(sample(20, 80, TRUE), sample(20, 80, TRUE)),
                             c(1, 1), c(1, 1), c(2, 3), c(2, 3)))
    p_flat <- rgraphtool:::cpp_sbm_check_move_prob(
      rgraphtool:::graph_pointer(g), sample(5, 20, TRUE), 0.5, 400L
    )
    expect_lt(max(p_flat), 1e-12)

    lv <- rgraphtool:::validate_hierarchy(
      list(sample(6, 20, TRUE), c(1, 1, 2, 2, 3, 3), c(1, 1, 2), c(1, 1)), 20
    )
    p_nested <- rgraphtool:::cpp_nested_check_move_prob(
      rgraphtool:::graph_pointer(g), lv, 1.0, 400L
    )
    expect_lt(max(p_nested), 1e-12)
  }
})

test_that("flat single-node and merge-split MCMC match exact Boltzmann distribution on N=3", {
  g <- graph_create(3, directed = FALSE)
  graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
  parts <- list("111" = c(1L, 1L, 1L), "112" = c(1L, 1L, 2L),
                "121" = c(1L, 2L, 1L), "122" = c(1L, 2L, 2L),
                "123" = c(1L, 2L, 3L))
  canon <- function(b) paste(match(b, unique(b)), collapse = "")
  freqs <- function(samples) {
    keys <- vapply(samples, canon, character(1))
    vapply(names(parts), function(nm) mean(keys == nm), numeric(1))
  }

  for (dc in c(FALSE, TRUE)) {
    S <- vapply(parts, function(b) block_description_length(block_state(g, b, dc)), numeric(1))
    exact <- exp(-S) / sum(exp(-S))
    st0 <- block_state(g, c(1L, 1L, 1L), degree_corrected = dc)

    set.seed(11)
    s_single <- sample_posterior(st0, nsamples = 25000L, thin = 1L, burnin = 300L,
                                 multiflip = FALSE, mcmc_args = list(beta = 1, c = 0.5, d = 0.2))
    expect_lt(max(abs(freqs(s_single$samples) - exact)), 0.02)

    set.seed(12)
    s_multi <- sample_posterior(st0, nsamples = 25000L, thin = 1L, burnin = 300L,
                                multiflip = TRUE,
                                mcmc_args = list(beta = 1, c = 0.5, d = 0.2, gibbs_sweeps = 3))
    expect_lt(max(abs(freqs(s_multi$samples) - exact)), 0.02)
  }

  # Pure split/merge/merge-split (psingle = 0) must also hit the exact stationary distribution
  S_ndc <- vapply(parts, function(b) block_description_length(block_state(g, b, FALSE)), numeric(1))
  exact_ndc <- exp(-S_ndc) / sum(exp(-S_ndc))
  st0 <- block_state(g, c(1L, 1L, 1L), degree_corrected = FALSE)
  set.seed(13)
  s_ms <- sample_posterior(
    st0, nsamples = 25000L, thin = 1L, burnin = 300L, multiflip = TRUE,
    mcmc_args = list(beta = 1, c = 0.5, d = 0.2, psingle = 0, psplit = 1,
                     pmerge = 1, pmergesplit = 1, gibbs_sweeps = 3)
  )
  expect_lt(max(abs(freqs(s_ms$samples) - exact_ndc)), 0.02)
})

test_that("directed flat MCMC matches exact Boltzmann distribution on N=3", {
  g <- graph_create(3, directed = TRUE)
  graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(1, 1)))
  parts <- list("111" = c(1L, 1L, 1L), "112" = c(1L, 1L, 2L),
                "121" = c(1L, 2L, 1L), "122" = c(1L, 2L, 2L),
                "123" = c(1L, 2L, 3L))
  canon <- function(b) paste(match(b, unique(b)), collapse = "")
  S <- vapply(parts, function(b) block_description_length(block_state(g, b, TRUE)), numeric(1))
  exact <- exp(-S) / sum(exp(-S))
  st0 <- block_state(g, c(1L, 1L, 1L), degree_corrected = TRUE)
  set.seed(14)
  s <- sample_posterior(st0, nsamples = 25000L, thin = 1L, burnin = 300L,
                        multiflip = TRUE, mcmc_args = list(beta = 1, c = 0.5, d = 0.2, gibbs_sweeps = 3))
  keys <- vapply(s$samples, canon, character(1))
  emp <- vapply(names(parts), function(nm) mean(keys == nm), numeric(1))
  expect_lt(max(abs(emp - exact)), 0.02)
})

test_that("nested merge-split MCMC matches exact stationary distribution over all 12 hierarchies", {
  g <- graph_create(3, directed = FALSE)
  graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
  all_parts <- function(n) {
    if (n == 1L) return(list(1L))
    if (n == 2L) return(list(c(1L, 1L), c(1L, 2L)))
    list(c(1L, 1L, 1L), c(1L, 1L, 2L), c(1L, 2L, 1L), c(1L, 2L, 2L), c(1L, 2L, 3L))
  }
  canon_h <- function(bs) paste(vapply(bs, function(b) paste(b, collapse = ""), character(1)), collapse = "/")
  S_h <- numeric()
  for (b0 in all_parts(3L)) {
    for (b1 in all_parts(max(b0))) {
      bs <- list(b0, b1)
      S_h[canon_h(bs)] <- nested_description_length(nested_block_state(g, bs, FALSE))
    }
  }
  exact <- exp(-S_h) / sum(exp(-S_h))
  nst0 <- nested_block_state(g, list(c(1L, 1L, 1L), 1L), degree_corrected = FALSE)

  set.seed(21)
  s0 <- sample_posterior(nst0, nsamples = 25000L, thin = 1L, burnin = 300L,
                         multiflip = TRUE, mcmc_args = list(beta = 0, c = 1, d = 0.3, gibbs_sweeps = 2))
  k0 <- vapply(s0$samples, canon_h, character(1))
  emp0 <- vapply(names(S_h), function(nm) mean(k0 == nm), numeric(1))
  expect_lt(max(abs(emp0 - 1 / 12)), 0.015)

  set.seed(22)
  s1 <- sample_posterior(nst0, nsamples = 25000L, thin = 1L, burnin = 300L,
                         multiflip = TRUE, mcmc_args = list(beta = 1, c = 1, d = 0.3, gibbs_sweeps = 2))
  k1 <- vapply(s1$samples, canon_h, character(1))
  emp1 <- vapply(names(S_h), function(nm) mean(k1 == nm), numeric(1))
  expect_lt(max(abs(emp1 - exact)), 0.015)
})

test_that("mcmc_sweep, mcmc_equilibrate, and mcmc_anneal respect bounds, seeds, and callbacks", {
  set.seed(30)
  g <- graph_create(30, directed = FALSE)
  edges <- rbind(
    do.call(rbind, lapply(list(1:10, 11:20, 21:30), function(v) t(combn(v, 2)))),
    c(1, 11), c(11, 21), c(21, 1)
  )
  graph_add_edges(g, edges)

  st <- block_state(g, rep(1:3, each = 10), degree_corrected = TRUE)
  set.seed(7)
  r1 <- mcmc_sweep(st, beta = 1, niter = 5L)
  set.seed(7)
  r2 <- mcmc_sweep(st, beta = 1, niter = 5L)
  expect_identical(block_membership(r1$state), block_membership(r2$state))
  expect_equal(r1$entropy, block_description_length(r1$state), tolerance = 1e-8)

  # Zero-temperature sweep never worsens description length
  st_rand <- block_state(g, sample(6, 30, TRUE), degree_corrected = TRUE)
  before <- block_description_length(st_rand)
  r_inf <- mcmc_sweep(st_rand, beta = Inf, niter = 5L, multiflip = TRUE)
  expect_lte(r_inf$entropy, before + 1e-8)

  # B_min / B_max bounds respected
  r_bnd <- mcmc_sweep(st, beta = 1, niter = 20L, B_min = 3L, B_max = 3L)
  expect_equal(length(block_sizes(r_bnd$state)), 3L)

  # Equilibration with history and collection, and R callback
  set.seed(8)
  eq <- mcmc_equilibrate(st, wait = 20L, nbreaks = 1L, max_niter = 40L,
                         history = TRUE, collect = TRUE, thin = 2L)
  expect_true(is.data.frame(eq$history))
  expect_equal(length(eq$samples), floor(eq$niter / 2))

  seen <- 0L
  eq_cb <- mcmc_equilibrate(st, force_niter = 5L, callback = function(s) {
    seen <<- seen + 1L
    length(block_sizes(s))
  })
  expect_equal(seen, 5L)
  expect_length(eq_cb$callback_values, 5L)

  # Nested annealing
  nst <- minimize_nested_blockmodel_dl(g)
  ann <- mcmc_anneal(nst, beta_range = c(1, 5), niter = 5L, history = TRUE,
                     mcmc_equilibrate_args = list(force_niter = 2L))
  expect_true(inherits(ann$state, "rgraphtool_nested_block_state"))
  expect_equal(ann$entropy, nested_description_length(ann$state), tolerance = 1e-8)
  expect_true(is.data.frame(ann$history))
  expect_equal(min(ann$history$beta), 1, tolerance = 1e-8)
  expect_equal(max(ann$history$beta), 5, tolerance = 1e-8)
})
