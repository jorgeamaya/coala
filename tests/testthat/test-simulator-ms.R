context('Simulator ms')

test_that("the ms sim program exists", {
  expect_that(get_simulator("ms"), is_a("Simulator"))
})


test_that("msSimFunc is working", {
  ms <- get_simulator("ms")
  dm_tt <- model_theta_tau()
  set.seed(789)
  sum_stats <- ms$simulate(dm_tt, c(1, 10))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
})


test_that("simulation works with recombination", {
  ms <- get_simulator("ms")
  dm_tt <- model_theta_tau() + feat_recombination(1)
  sum_stats <- ms$simulate(dm_tt, c(1, 10))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
})


test_that("msSimFunc works with inter-locus variation", {
  ms <- get_simulator("ms")
  dm_tmp <- coal_model(5:6, 2) +
    feat_mutation(par_const(2), variance = '2') +
    feat_pop_merge(par_const(.5), 2, 1) +
    sumstat_jsfs()

  set.seed(117)
  sum_stats <- ms$simulate(dm_tmp)
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)

  set.seed(117)
  sum_stats2 <- ms$simulate(dm_tmp)
  expect_equal(sum_stats$jsfs, sum_stats2$jsfs)
})


test_that("simulating a size change works", {
  ms <- get_simulator("ms")
  dm_tmp <- coal_model(5:6, 1) +
    feat_mutation(par_range('theta', 1, 10), variance = '2') +
    feat_pop_merge(par_const(.5), 2, 1) +
    feat_size_change(par_const(2), population = 2, at.time = 0) +
    sumstat_jsfs()

  sum_stats <- ms$simulate(dm_tmp, 5)
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
})


test_that('Printing the command works', {
  ms <- get_simulator("ms")
  cmd <- ms$get_cmd(coal_model(5:6, 17))
  expect_equal(grep('ms', cmd), 1)
  expect_equal(grep('ms 11 17', cmd), 1)
  expect_equal(grep('-I 2 5 6', cmd), 1)

  cmd <- ms$get_cmd(coal_model(5:6, 17) +
                          par_range('a', 1, 5) +
                          feat_mutation(par_expr(2 * a)))
  expect_equal(grep('ms', cmd), 1)
  expect_equal(grep('ms 11 17', cmd), 1)
  expect_equal(grep('-I 2 5 6', cmd), 1)
  expect_equal(grep('a', cmd), 1)
})


test_that('simulating unphased data works', {
  ms <- get_simulator("ms")
  model <- model_theta_tau() + feat_unphased(2, 1) + sumstat_seg_sites()
  stats <- ms$simulate(model, c(1,5))
  expect_equal(dim(stats$jsfs), c(11, 16))
  expect_equal(nrow(stats$seg_sites[[1]]), 25)

  model <- model_theta_tau() + feat_unphased(3, 2) + sumstat_seg_sites()
  stats <- ms$simulate(model, c(1,5))
  expect_equal(dim(stats$jsfs), c(21, 31))
  expect_equal(nrow(stats$seg_sites[[1]]), 50)
})