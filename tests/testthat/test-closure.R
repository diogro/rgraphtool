# Tests for Milestone 10: LatentClosureBlockState (Peixoto 2022 PRX)

make_closure_ref_graph <- function() {
  g <- graph_create(5L, directed = FALSE)
  e5 <- matrix(c(1L,2L, 2L,3L, 3L,4L, 4L,5L, 1L,3L, 3L,5L, 1L,5L), ncol = 2L, byrow = TRUE)
  graph_add_edges(g, e5)
  ge <- graph_edges(g)
  ek <- paste0(pmin(ge[, 1L], ge[, 2L]), ":", pmax(ge[, 1L], ge[, 2L]))
  gen <- ifelse(ek %in% c("1:3", "3:5"), 1L, ifelse(ek == "1:5", 2L, 0L))
  ex <- lapply(ek, function(k) {
    if (k == "1:3") 2L else if (k == "3:5") 4L else if (k == "1:5") 3L else integer(0)
  })
  list(g = g, gen = gen, ex = ex)
}

test_that("latent_closure_block_state matches pinned graph-tool 3.9 for L=2 flat & nested", {
  obj <- make_closure_ref_graph()
  b5 <- c(1L, 1L, 2L, 2L, 2L)
  bs5 <- list(b5, c(1L, 1L))

  # Flat L=2
  st_f <- latent_closure_block_state(
    obj$g, L = 2L, nested = FALSE, blocks = b5, gen = obj$gen, ex = obj$ex
  )
  dl_f <- closure_description_length(st_f, breakdown = TRUE)
  expect_equal(dl_f$total, 19.203211485647, tolerance = 1e-11)
  expect_equal(dl_f$sbm, 13.827933077962832, tolerance = 1e-11)
  expect_equal(dl_f$layer_entropies[1L], 2.484906649788, tolerance = 1e-11)
  expect_equal(dl_f$layer_entropies[2L], 2.8903717578961645, tolerance = 1e-11)
  expect_equal(dl_f$M[[1L]], c(0, 1, 1, 1, 0))
  expect_equal(dl_f$E[[1L]], c(0, 1, 0, 1, 0))
  expect_equal(dl_f$M[[2L]], c(0, 0, 3, 0, 0))
  expect_equal(dl_f$E[[2L]], c(0, 0, 1, 0, 0))
  expect_equal(closure_description_length(st_f), 19.203211485647, tolerance = 1e-11)
  expect_output(print(st_f), "rgraphtool_latent_closure_block_state")

  # Nested L=2
  st_n <- latent_closure_block_state(
    obj$g, L = 2L, nested = TRUE, blocks = bs5, gen = obj$gen, ex = obj$ex
  )
  dl_n <- closure_description_length(st_n, breakdown = TRUE)
  expect_equal(dl_n$total, 19.89635866620694, tolerance = 1e-11)
  expect_equal(dl_n$sbm, 14.521080258522776, tolerance = 1e-11)
  expect_equal(dl_n$layer_entropies[1L], 2.484906649788, tolerance = 1e-11)
  expect_equal(dl_n$layer_entropies[2L], 2.8903717578961645, tolerance = 1e-11)
})

test_that("multi-ego triadic closure matches exact entropy", {
  # Diamond 1-2, 1-3, 2-4, 3-4 plus chord 2-3 closed by both egos 1 and 4
  g <- graph_create(4L, directed = FALSE)
  graph_add_edges(g, matrix(c(1L,2L, 1L,3L, 2L,4L, 3L,4L, 2L,3L), ncol = 2L, byrow = TRUE))
  ge <- graph_edges(g)
  ek <- paste0(pmin(ge[, 1L], ge[, 2L]), ":", pmax(ge[, 1L], ge[, 2L]))
  gen <- ifelse(ek == "2:3", 1L, 0L)
  ex <- lapply(ek, function(k) if (k == "2:3") c(1L, 4L) else integer(0))
  st <- latent_closure_block_state(g, L = 1L, nested = FALSE, blocks = rep(1L, 4L),
                                   gen = gen, ex = ex)
  dl <- closure_description_length(st, breakdown = TRUE)
  expect_equal(dl$M[[1L]], c(1, 1, 1, 1))
  expect_equal(dl$E[[1L]], c(1, 0, 0, 1))
  # N_M = 4, N_E = 2 -> log(5) + lchoose(4, 2) + 0 + 0 = log(30)
  expect_equal(dl$closure, log(30), tolerance = 1e-12)
})

test_that("closure incremental moves match full recomputation", {
  set.seed(505)
  obj <- make_closure_ref_graph()
  b5 <- c(1L, 1L, 2L, 2L, 2L)
  bs5 <- list(b5, c(1L, 1L))
  for (nest in c(FALSE, TRUE)) {
    st <- latent_closure_block_state(
      obj$g, L = 2L, nested = nest,
      blocks = if (nest) bs5 else b5, gen = obj$gen, ex = obj$ex
    )
    errs <- rgraphtool:::cpp_closure_check_moves(
      st$n_vertices, st$L, st$obs_edges, st$gen, st$w0, st$loops, st$ex,
      st$bs, st$degree_corrected, st$dist, 120L, Inf
    )
    expect_lt(max(errs), 1e-11)
  }
})

test_that("infer_triadic_closure distinguishes triangle-closing edges from substrate edges", {
  set.seed(777)
  # Sparse 20-cycle + 2 long chords (triangle-free substrate) + 10 triadic closures
  g <- graph_create(20L, directed = FALSE)
  sub_e <- rbind(cbind(1:19, 2:20), c(20L, 1L), c(1L, 10L), c(5L, 15L))
  cl_e <- cbind(seq(1L, 19L, by = 2L), seq(3L, 21L, by = 2L) %% 20L)
  cl_e[cl_e == 0L] <- 20L
  graph_add_edges(g, rbind(sub_e, cl_e))

  st <- latent_closure_block_state(g, L = 1L, nested = FALSE)
  res <- infer_triadic_closure(st, nsamples = 200L, burnin = 150L, thin = 2L,
                               init_sbm = FALSE, collect_partitions = TRUE)
  expect_s3_class(res, "rgraphtool_closure_inference")
  expect_length(res$partitions, 200L)
  expect_equal(nrow(res$edges), graph_num_edges(g))
  expect_equal(res$edges$substrate_prob + res$edges$closure_prob,
               rep(1, nrow(res$edges)), tolerance = 1e-12)
  expect_output(print(res), "rgraphtool_closure_inference")

  ek <- function(u, v) paste0(pmin(u, v), ":", pmax(u, v))
  is_cl <- ek(res$edges$from, res$edges$to) %in% ek(cl_e[, 1L], cl_e[, 2L])
  expect_gt(mean(res$edges$closure_prob[is_cl]), 0.65)
  expect_lt(mean(res$edges$closure_prob[!is_cl]), 0.25)

  # Even vertices are the closing egos
  vs1 <- res$vertex_stats$gen_1
  expect_gt(mean(vs1$t[seq(2L, 20L, by = 2L)]), mean(vs1$t[seq(1L, 19L, by = 2L)]))
})

test_that("latent_closure_block_state validates invalid inputs and states", {
  g_dir <- graph_create(4L, directed = TRUE)
  graph_add_edges(g_dir, matrix(c(1L, 2L), ncol = 2L))
  expect_error(latent_closure_block_state(g_dir), "undirected")

  g_loop <- graph_create(3L, directed = FALSE)
  graph_add_edges(g_loop, matrix(c(1L, 1L, 1L, 2L), ncol = 2L, byrow = TRUE))
  expect_error(latent_closure_block_state(g_loop), "self-loops")

  g_par <- graph_create(3L, directed = FALSE)
  graph_add_edges(g_par, matrix(c(1L, 2L, 1L, 2L), ncol = 2L, byrow = TRUE))
  expect_error(latent_closure_block_state(g_par), "parallel")

  # Invalid ego (not forming a triangle)
  g_path <- graph_create(3L, directed = FALSE)
  graph_add_edges(g_path, matrix(c(1L, 2L, 2L, 3L), ncol = 2L, byrow = TRUE))
  expect_error(
    latent_closure_block_state(g_path, L = 1L, gen = c(0L, 1L), ex = list(integer(0), 1L)),
    "Invalid initial state"
  )
})
