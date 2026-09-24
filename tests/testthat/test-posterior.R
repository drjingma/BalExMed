make_problem <- function(seed = 1, n = 40, d = 6) {
  set.seed(seed)
  sim <- simulate_balexmed(n = n, d = d, n_plus = 2, n_minus = 1, alpha = 0.8,
                           beta = 0.6, gamma = 0.7)
  C <- cbind(1, sim$covariates$c, rnorm(n))
  list(X = sim$X, m = sim$m, y = sim$y, C = C)
}

test_that("fast log posterior equals the naive computation", {
  p <- make_problem()
  h <- 0.3; eta <- 0.2; nu <- 40; lm <- var(p$m); ly <- var(p$y)
  s <- BalExMed:::.bm_setup(p$X, p$m, p$y, p$C, eta, h, nu, lm, ly)
  zs <- list(c(1, 1, -1, 0, 0, 0), c(1, 0, -1, -1, 0, 1), c(-1, 1, 0, 0, 1, -1))
  for (z in zs) {
    B <- balance(p$X, z)
    fast <- BalExMed:::.bm_logpost(B, sum(z == 1), sum(z == -1), s)
    ref <- naive_logpost(z, p$X, p$m, p$y, p$C, eta, h, nu, lm, ly)
    expect_equal(fast, ref, tolerance = 1e-9)
  }
})

test_that("log posterior differences equal those of the exact marginal density", {
  p <- make_problem(seed = 4)
  h <- 0.5; eta <- 1/3; nu <- 40; lm <- var(p$m); ly <- var(p$y)
  s <- BalExMed:::.bm_setup(p$X, p$m, p$y, p$C, eta, h, nu, lm, ly)
  zs <- all_valid_z(6)[c(1, 50, 200, 400), ]
  fast <- apply(zs, 1, function(z)
    BalExMed:::.bm_logpost(balance(p$X, z), sum(z == 1), sum(z == -1), s))
  exact <- apply(zs, 1, function(z)
    exact_logpost(z, p$X, p$m, p$y, p$C, eta, h, nu, lm, ly))
  expect_equal(fast - fast[1], exact - exact[1], tolerance = 1e-8)
})

test_that("log posterior is symmetric under exchanging numerator and denominator", {
  p <- make_problem(seed = 5)
  s <- BalExMed:::.bm_setup(p$X, p$m, p$y, p$C, 0.2, 1e-6, 40, var(p$m), var(p$y))
  z <- c(1, 1, -1, 0, -1, 0)
  lp <- function(z) BalExMed:::.bm_logpost(balance(p$X, z), sum(z == 1), sum(z == -1), s)
  expect_equal(lp(z), lp(-z), tolerance = 1e-10)
})
