# Tests for Milestone 11: Network reconstruction from dynamics and behavior,
# L1 regularization, quantized MDL priors, and native forward simulators.

make_m11_graph <- function(directed, w) {
  g <- graph_create(4L, directed = directed)
  graph_add_edges(g, matrix(c(1L,2L, 2L,3L, 3L,4L), ncol = 2L, byrow = TRUE),
                  properties = list(w = w))
  g
}

test_that("all 10 dynamics models and priors match pinned graph-tool 3.9 to machine precision", {
  s_ising <- matrix(c(
    1,-1, 1, 1,-1,
   -1, 1, 1,-1, 1,
    1, 1,-1, 1,-1,
   -1,-1, 1,-1, 1
  ), nrow = 4L, byrow = TRUE)
  s_iz <- s_ising; s_iz[1, 2] <- 0; s_iz[3, 4] <- 0

  # Exact same s_norm from numpy seed 42 used in graph-tool 3.9 fixture
  s_norm <- matrix(c(
     0.4967141530112327, -0.13826430117118466,  0.6476885381006925,  1.5230298564080254, -0.23415337472333597,
    -0.23413695694918055, 1.5792128155073915,   0.7674347291529088, -0.4694743859349521,  0.5425600435859647,
    -0.46341769281246226,-0.46572975357025687,  0.24196227156603412,-1.913280244657798,  -1.7249178325130328,
    -0.5622875292409727, -1.0128311203344238,   0.3142473325952739, -0.9080240755212109, -1.4123037013352915,
     1.465648768921554,  -0.22577630048653566,  0.06752820468792384,-1.4247481862134568, -0.5443827245251827,
     0.11092258970986608,-1.1509935774223028,   0.37569801834567196,-0.600638689918805,  -0.2916937497932768,
    -0.6017066122293969,  1.8522781845089378,  -0.013497224737933921,-1.0577109289559004, 0.822544912103189,
    -1.2208436499710222,  0.2088635950047554,  -1.9596701238797756, -1.3281860488984305,  0.19686123586912352
  ), nrow = 4L, byrow = TRUE)
  s_ci <- tanh(s_norm[, 1:5])

  s_epi <- matrix(c(
    1,1,2,2,2,
    0,1,1,2,2,
    0,0,1,1,2,
    0,0,0,1,1
  ), nrow = 4L, byrow = TRUE)
  s_potts <- matrix(c(
    0,1,2,0,1,2,
    1,2,0,1,2,0,
    2,0,1,2,0,1,
    0,2,1,0,2,1
  ), nrow = 4L, byrow = TRUE)

  th1   <- matrix(c(0.1, -0.2, 0.3, 0.0), ncol = 1L)
  th_e  <- matrix(c(-2.0, -1.5, -2.5, -1.0), ncol = 1L)
  th_pn <- matrix(c(-0.2, 0.1, 0.1, -0.3, -0.1, 0.2, 0.0, -0.1), nrow = 4L, byrow = TRUE)
  th_ln <- matrix(c(-0.2, 0.1, -0.1, 0.0), ncol = 1L)
  th_pt <- matrix(c(0.1,-0.1,0.2, 0.0,0.2,-0.2, -0.1,0.1,0.0, 0.2,0.0,-0.1), nrow = 4L, byrow = TRUE)
  f_pp  <- pmax(matrix(c(1.0, 0.2, -0.1, 0.2, 1.0, 0.3, -0.1, 0.3, 1.0), nrow = 3L, byrow = TRUE), 0) / 3.5
  b4    <- c(1L, 1L, 2L, 2L)

  g_u   <- make_m11_graph(FALSE, c(0.5, -0.3, 0.8))
  g_d   <- make_m11_graph(TRUE,  c(0.5, -0.3, 0.8))
  g_e   <- make_m11_graph(TRUE,  c(-0.5, -0.8, -0.3))
  g_nu  <- make_m11_graph(FALSE, c(0.2, -0.15, 0.25))
  g_nup <- make_m11_graph(FALSE, c(1.5, -1.2, 1.8))
  g_nd  <- make_m11_graph(TRUE,  c(0.2, -0.15, 0.25))

  st_pi  <- pseudo_ising_block_state(s_ising, graph = g_u, x = "w", theta = th1, xdelta = 0.01, tdelta = 0.01, nested = FALSE, blocks = b4)
  st_piz <- pseudo_ising_block_state(s_iz, graph = g_u, has_zero = TRUE, x = "w", theta = th1, xdelta = 0.01, tdelta = 0.01, nested = FALSE, blocks = b4)
  st_ig  <- ising_glauber_block_state(s_ising, graph = g_d, x = "w", theta = th1, xdelta = 0.01, tdelta = 0.01, nested = FALSE, blocks = b4)
  st_pci <- pseudo_cising_block_state(s_ci, graph = g_u, x = "w", theta = th1)
  st_cig <- cising_glauber_block_state(s_ci, graph = g_d, x = "w", theta = th1)
  st_ep  <- epidemics_block_state(s_epi, graph = g_e, x = "w", theta = th_e)
  st_pn  <- pseudo_normal_block_state(s_norm, graph = g_nu, positive = FALSE, x = "w", theta = th_pn)
  st_pnp <- pseudo_normal_block_state(s_norm, graph = g_nup, positive = TRUE, pslack = 1e-6, x = "w", theta = th_pn)
  st_ng  <- normal_glauber_block_state(s_norm, graph = g_nd, x = "w", theta = th_pn)
  st_ln  <- linear_normal_block_state(s_norm, graph = g_nd, x = "w", theta = th_ln)
  st_pp  <- pseudo_potts_block_state(s_potts, f = f_pp, graph = g_u, x = "w", theta = th_pt)

  expect_equal(dynamics_description_length(st_pi,  breakdown = TRUE)$nll, 26.576890671264938, tolerance = 1e-11)
  expect_equal(dynamics_description_length(st_piz, breakdown = TRUE)$nll, 31.354818125662714, tolerance = 1e-11)
  expect_equal(dynamics_description_length(st_ig,  breakdown = TRUE)$nll, 11.773898998251784, tolerance = 1e-11)
  expect_equal(dynamics_description_length(st_pci, breakdown = TRUE)$nll, 15.971484693063514, tolerance = 1e-11)
  expect_equal(dynamics_description_length(st_cig, breakdown = TRUE)$nll, 12.106584938835596, tolerance = 1e-11)
  expect_equal(dynamics_description_length(st_ep,  breakdown = TRUE)$nll,  2.803574130345813, tolerance = 1e-11)
  expect_equal(dynamics_description_length(st_pn,  breakdown = TRUE)$nll, 51.082291667536836, tolerance = 1e-11)
  expect_equal(dynamics_description_length(st_pnp, breakdown = TRUE)$nll, 63.945008827502143, tolerance = 1e-11)
  expect_equal(dynamics_description_length(st_ng,  breakdown = TRUE)$nll, 50.103361537558385, tolerance = 1e-11)
  expect_equal(dynamics_description_length(st_ln,  breakdown = TRUE)$nll, 63.397511540161901, tolerance = 1e-11)
  expect_equal(dynamics_description_length(st_pp,  breakdown = TRUE)$nll, 27.134386666017619, tolerance = 1e-11)

  # Exact priors matching graph-tool 3.9
  bd_sbm <- dynamics_description_length(st_pi, sbm = TRUE,  xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 1.0, breakdown = TRUE)
  bd_l1  <- dynamics_description_length(st_pi, sbm = FALSE, xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 1.0, breakdown = TRUE)
  bd_qd  <- dynamics_description_length(st_pi, sbm = FALSE, xdist = TRUE,  tdist = TRUE,  xl1 = 1.0, tl1 = 1.0, breakdown = TRUE)
  bd_dir <- dynamics_description_length(st_ig, sbm = FALSE, xdist = FALSE, tdist = FALSE, xl1 = 0.0, tl1 = 0.0, breakdown = TRUE)

  expect_equal(bd_sbm$geom_prior, 2.408362188537088,  tolerance = 1e-12)
  expect_equal(bd_dir$geom_prior, 2.805077480482146,  tolerance = 1e-12)
  expect_equal(bd_l1$x_prior,    33.397388574299306,  tolerance = 1e-11)
  expect_equal(bd_l1$t_prior,    21.795755924536131,  tolerance = 1e-11)
  expect_equal(bd_qd$x_prior,   146.565982204230181,  tolerance = 1e-11)
  expect_equal(bd_qd$t_prior,   150.037701539136179,  tolerance = 1e-11)
  expect_equal(bd_sbm$total,     78.423718534458914,  tolerance = 1e-11)
  expect_output(print(st_pi), "rgraphtool_dynamics_block_state")
})

test_that("incremental edge, node, and SBM moves match full recomputation across all models", {
  set.seed(1101)
  s_ising <- matrix(sample(c(-1L, 1L), 24L, replace = TRUE), nrow = 4L)
  s_cont  <- matrix(stats::rnorm(32L), nrow = 4L)
  s_ci    <- tanh(s_cont[, 1:6])
  s_epi   <- matrix(c(1,1,2,2,2, 0,1,1,2,2, 0,0,1,1,2, 0,0,0,1,1), nrow = 4L, byrow = TRUE)
  s_potts <- matrix(sample(0:2, 24L, replace = TRUE), nrow = 4L)

  g_u <- make_m11_graph(FALSE, c(0.4, -0.3, 0.5))
  g_d <- make_m11_graph(TRUE,  c(0.4, -0.3, 0.5))
  g_e <- make_m11_graph(TRUE,  c(-0.4, -0.6, -0.3))
  b4  <- c(1L, 1L, 2L, 2L)

  states <- list(
    pseudo_ising_block_state(s_ising, graph = g_u, x = "w", blocks = b4),
    ising_glauber_block_state(s_ising, graph = g_d, x = "w", blocks = b4),
    pseudo_cising_block_state(s_ci, graph = g_u, x = "w", blocks = b4),
    cising_glauber_block_state(s_ci, graph = g_d, x = "w", blocks = b4),
    epidemics_block_state(s_epi, graph = g_e, x = "w", blocks = b4),
    pseudo_normal_block_state(s_cont, graph = g_u, positive = TRUE, x = "w", blocks = b4),
    normal_glauber_block_state(s_cont, graph = g_d, x = "w", blocks = b4),
    linear_normal_block_state(s_cont, graph = g_d, x = "w", blocks = b4),
    pseudo_potts_block_state(s_potts, f = diag(3), graph = g_u, x = "w", blocks = b4),
    potts_glauber_block_state(s_potts, f = diag(3), graph = g_d, x = "w", blocks = b4)
  )

  for (st in states) {
    for (qd in c(FALSE, TRUE)) {
      opts <- st$default_prior
      opts$sbm <- TRUE; opts$xdist <- qd; opts$tdist <- rep(qd, st$q_t)
      opts$xl1 <- 1.0; opts$tl1 <- rep(1.0, st$q_t)
      opts$xdelta <- 0.01; opts$tdelta <- rep(0.01, st$q_t)
      errs <- rgraphtool:::cpp_dynamics_check_moves(
        st$model_code, st$n_vertices, st$directed, st$self_loops,
        st$s_list, st$f_mat, st$has_zero, st$exposed,
        st$positive, st$pslack, st$x_range[1L], st$x_range[2L],
        st$tmin, st$tmax, st$edges, st$x, st$theta, st$bs,
        st$degree_corrected, st$dist, opts, 30L, Inf
      )
      expect_lt(max(errs), 1e-10)
    }
  }
})

test_that("Ising, Epidemics, Normal, Linear, and Potts reconstruction recover true graphs", {
  set.seed(2026)
  ek <- function(u, v, dir = FALSE) if (!dir) paste0(pmin(u, v), ":", pmax(u, v)) else paste0(u, ":", v)

  # 1. Undirected sparse 2-ring graph (N = 12, 13 edges)
  g_u <- graph_create(12L, directed = FALSE)
  e_u <- rbind(cbind(1:5, 2:6), c(6L, 1L), cbind(7:11, 8:12), c(12L, 7L), c(1L, 7L))
  graph_add_edges(g_u, e_u)
  ap_u <- t(utils::combn(12L, 2L))
  y_u <- as.integer(ek(ap_u[, 1L], ap_u[, 2L]) %in% ek(e_u[, 1L], e_u[, 2L]))

  # 1a. PseudoIsing (SBM + L1 + posterior MCMC)
  s_pi <- simulate_ising(g_u, weights = 0.6, n_samples = 450L, burnin = 200L, thin = 4L, mode = "equilibrium")
  st_pi <- pseudo_ising_block_state(s_pi, nested = TRUE, xdelta = 1e-4, tdelta = 1e-4)
  rec_pi <- reconstruct_dynamics(st_pi, niter = 15L, nsamples = 20L, burnin = 10L, beta = 1,
                                 sbm = TRUE, xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 0.1)
  expect_s3_class(rec_pi, "rgraphtool_dynamics_reconstruction")
  expect_output(print(rec_pi), "rgraphtool_dynamics_reconstruction")
  idx <- match(ek(ap_u[, 1L], ap_u[, 2L]), ek(rec_pi$edges$from, rec_pi$edges$to))
  sc <- ifelse(is.na(idx), 0, abs(rec_pi$edges$weight[idx]) * rec_pi$edges$prob[idx])
  expect_gt(evaluate_edge_prediction(sc, y_u)$roc_auc, 0.95)

  # 1b. L1 regularization path
  path <- reconstruct_l1_path(st_pi, xl1_values = c(1, 5, 15), sbm = FALSE, niter = 10L, truth = g_u)
  expect_s3_class(path, "rgraphtool_l1_path")
  expect_equal(nrow(path$summary), 3L)
  expect_true(path$summary$n_edges[1L] <= path$summary$n_edges[3L] + 2L)
  expect_true(all(path$summary$roc_auc > 0.90))
  expect_output(print(path), "rgraphtool_l1_path")

  # 2. Epidemics multi-cascade SIS reconstruction (directed)
  g_d <- graph_create(10L, directed = TRUE)
  e_d <- rbind(cbind(1:9, 2:10), c(10L, 1L), c(1L, 6L), c(6L, 1L))
  graph_add_edges(g_d, e_d)
  casc <- simulate_epidemics(g_d, transmission = 0.45, spontaneous = 0.01, recovery = 0.3,
                             model = "SIS", n_steps = 50L, n_cascades = 15L)
  st_ep <- epidemics_block_state(casc, directed = TRUE, nested = FALSE, xdelta = 1e-4, tdelta = 1e-4)
  rec_ep <- reconstruct_dynamics(st_ep, niter = 12L, sbm = TRUE, xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 0.1)
  tp_ep <- dynamics_transmission_probs(rec_ep)
  expect_true(nrow(tp_ep$edges) >= 10L)
  expect_true(all(tp_ep$edges$transmission_prob > 0 & tp_ep$edges$transmission_prob < 1))
  ap_d <- expand.grid(u = 1:10, v = 1:10); ap_d <- ap_d[ap_d$u != ap_d$v, ]
  y_d <- as.integer(ek(ap_d$u, ap_d$v, TRUE) %in% ek(e_d[, 1L], e_d[, 2L], TRUE))
  idx_d <- match(ek(ap_d$u, ap_d$v, TRUE), ek(rec_ep$edges$from, rec_ep$edges$to, TRUE))
  sc_d <- ifelse(is.na(idx_d), 0, abs(rec_ep$edges$weight[idx_d]))
  expect_gt(evaluate_edge_prediction(sc_d, y_d)$roc_auc, 0.95)

  # 3. PseudoNormal precision matrix reconstruction
  s_pn <- simulate_normal_dynamics(g_u, weights = 0.85, sigma = 0.6, mu = 0,
                                   n_samples = 600L, burnin = 200L, thin = 3L, mode = "equilibrium")
  st_pn <- pseudo_normal_block_state(s_pn, fix_mean = TRUE, positive = TRUE, xdelta = 1e-4, tdelta = 1e-4)
  rec_pn <- reconstruct_dynamics(st_pn, niter = 15L, sbm = TRUE, xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 0.1)
  prec <- dynamics_precision_matrix(rec_pn)
  expect_true(min(eigen(prec$precision, symmetric = TRUE, only.values = TRUE)$values) > 0)
  idx_pn <- match(ek(ap_u[, 1L], ap_u[, 2L]), ek(rec_pn$edges$from, rec_pn$edges$to))
  sc_pn <- ifelse(is.na(idx_pn), 0, abs(rec_pn$edges$weight[idx_pn]))
  expect_gt(evaluate_edge_prediction(sc_pn, y_u)$roc_auc, 0.95)

  # 4. Potts categorical reconstruction
  s_pt <- simulate_potts(g_u, weights = 1.2, f = diag(3), n_samples = 350L, burnin = 150L, thin = 4L, mode = "equilibrium")
  st_pt <- pseudo_potts_block_state(s_pt, f = diag(3), nested = FALSE, xdelta = 1e-4, tdelta = 1e-4)
  rec_pt <- reconstruct_dynamics(st_pt, niter = 12L, sbm = TRUE, xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 0.1)
  idx_pt <- match(ek(ap_u[, 1L], ap_u[, 2L]), ek(rec_pt$edges$from, rec_pt$edges$to))
  sc_pt <- ifelse(is.na(idx_pt), 0, abs(rec_pt$edges$weight[idx_pt]))
  expect_gt(evaluate_edge_prediction(sc_pt, y_u)$roc_auc, 0.95)
})

test_that("dynamics constructors and helpers validate inputs cleanly", {
  s_ok <- matrix(c(1, -1, -1, 1, 1, 1), nrow = 3L)
  expect_error(pseudo_ising_block_state("bad"), "matrix")
  expect_error(pseudo_ising_block_state(matrix(c(1, NA, 1, 1), nrow = 2L)), "finite")
  expect_error(pseudo_ising_block_state(s_ok, x_range = c(2, 1)), "x_range")
  expect_error(pseudo_potts_block_state(matrix(c(0, 1, 5, 0), nrow = 2L), q = 3L), "0\\.\\.2")
  expect_error(dynamics_precision_matrix(pseudo_ising_block_state(s_ok)), "pseudo_normal")
  expect_error(dynamics_transmission_probs(pseudo_ising_block_state(s_ok)), "epidemics")

  g <- graph_create(3L, directed = FALSE)
  graph_add_edges(g, matrix(c(1L, 2L), ncol = 2L))
  expect_error(simulate_epidemics(g, transmission = 1.2), "\\[0, 1\\)")
  expect_error(simulate_normal_dynamics(g, sigma = -1), "positive")
})
