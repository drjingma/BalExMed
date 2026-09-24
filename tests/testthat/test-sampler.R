exact_posterior <- function(X, m, y, covariates, eta, h = 1e-6) {
  n <- nrow(X)
  C <- cbind(1, as.matrix(covariates))
  nu <- n; lm <- var(m); ly <- var(y)
  Z <- all_valid_z(ncol(X))
  lp <- apply(Z, 1, function(z) naive_logpost(z, X, m, y, C, eta, h, nu, lm, ly))
  post <- exp(lp - max(lp)); post <- post / sum(post)
  g <- apply(Z, 1, function(z) prob_gamma_positive(z, X, m, y, C, h, nu, ly))
  post_relabeled <- 2 * post * g          # distribution of z under gamma > 0
  list(PIP = colSums(post * (Z != 0)),
       P_plus = colSums(post_relabeled * (Z == 1)))
}

test_that("Gibbs and Metropolis-Hastings recover the exact posterior of z", {
  skip_on_cran()
  ## weak signal, so that several inclusion probabilities are far from 0 and 1
  set.seed(3)
  sim <- simulate_balexmed(n = 30, d = 5, n_plus = 1, n_minus = 1, alpha = 0.3,
                           beta = 0.3, gamma = 0.25, var_log = 1)
  ex <- exact_posterior(sim$X, sim$m, sim$y, sim$covariates, eta = 1/3)
  expect_gt(sum(ex$PIP > 0.15 & ex$PIP < 0.85), 2)
  fit_g <- balexmed(sim$X, sim$m, sim$y, sim$covariates, eta = 1/3,
                    n_iter = 20000, burn_in = 1000, seed = 21, verbose = FALSE)
  fit_mh <- balexmed(sim$X, sim$m, sim$y, sim$covariates, eta = 1/3, sampler = "mh",
                     n_iter = 150000, burn_in = 1000, seed = 21, verbose = FALSE)
  for (fit in list(fit_g, fit_mh)) {
    tt <- taxa_table(fit)
    tt <- tt[match(colnames(sim$X), tt$taxon), ]
    expect_lt(max(abs(tt$PIP - ex$PIP)), 0.03)
    expect_lt(max(abs(tt$P_plus - ex$P_plus)), 0.03)
  }
})

test_that("draws satisfy the sign convention and seeds give reproducible chains", {
  set.seed(7)
  sim <- simulate_balexmed(n = 80, d = 12)
  fit1 <- balexmed(sim$X, sim$m, sim$y, sim$covariates, n_iter = 300, burn_in = 100,
                   n_chains = 2, seed = 5, verbose = FALSE)
  fit2 <- balexmed(sim$X, sim$m, sim$y, sim$covariates, n_iter = 300, burn_in = 100,
                   n_chains = 2, seed = 5, verbose = FALSE)
  expect_true(all(effect_draws(fit1)$gamma > 0))
  expect_identical(fit1$chains[[2]]$z, fit2$chains[[2]]$z)
  expect_identical(fit1$chains[[1]]$alpha, fit2$chains[[1]]$alpha)
  skip_on_os("windows")
  fit3 <- balexmed(sim$X, sim$m, sim$y, sim$covariates, n_iter = 300, burn_in = 100,
                   n_chains = 2, seed = 5, cores = 2, verbose = FALSE)
  expect_identical(fit1$chains[[2]]$gamma, fit3$chains[[2]]$gamma)
})

test_that("the sampler recovers a simulated balance and its effects", {
  skip_on_cran()
  set.seed(8)
  sim <- simulate_balexmed(n = 200, d = 20)
  fit <- balexmed(sim$X, sim$m, sim$y, sim$covariates, n_iter = 3000, burn_in = 1000,
                  seed = 3, verbose = FALSE)
  expect_equal(unname(selected_z(fit)), unname(sim$z))
  sm <- summary(fit)$effects
  expect_true(sm["gamma (direct effect)", "2.5%"] < 1/3 &&
                sm["gamma (direct effect)", "97.5%"] > 1/3)
  expect_true(sm["alpha*beta (indirect effect)", "2.5%"] < 2/3 &&
                sm["alpha*beta (indirect effect)", "97.5%"] > 2/3)
})
