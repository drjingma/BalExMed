test_that("balexmed returns a well-formed fit and its methods work", {
  set.seed(9)
  sim <- simulate_balexmed(n = 60, d = 10)
  cov <- data.frame(c = sim$covariates$c, group = factor(rep(c("a", "b", "c"), 20)))
  fit <- balexmed(sim$X, sim$m, sim$y, covariates = cov, n_iter = 200, burn_in = 100,
                  thin = 2, n_chains = 2, z_init = "pca", seed = 1, verbose = FALSE)
  expect_s3_class(fit, "balexmed")
  expect_equal(fit$k, 4)
  expect_equal(dim(fit$chains[[1]]$z), c(50, 10))
  expect_equal(colnames(fit$chains[[1]]$psi_m), c("(Intercept)", "c", "groupb", "groupc"))
  expect_equal(unname(fit$chains[[1]]$z_init), pca_start(sim$X))

  dr <- effect_draws(fit)
  expect_equal(nrow(dr), 100)
  expect_equal(dr$indirect, dr$alpha * dr$beta)

  tt <- taxa_table(fit)
  expect_equal(nrow(tt), 10)
  expect_true(all(abs(tt$P_plus + tt$P_minus - tt$PIP) < 1e-12))
  expect_true("PIP_sd_chains" %in% names(tt))

  sm <- summary(fit)
  expect_equal(dim(sm$effects), c(5, 7))
  expect_output(print(sm), "indirect effect")
  expect_output(print(fit), "Balance mediation fit")
  expect_s3_class(as_mcmc_list(fit), "mcmc.list")

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(fit))
  expect_invisible(plot(fit, type = "inclusion", top = 5))
})

test_that("covariate matrices gain an intercept only when they lack one", {
  set.seed(10)
  sim <- simulate_balexmed(n = 50, d = 8)
  C1 <- cbind(x = rnorm(50))
  fit1 <- balexmed(sim$X, sim$m, sim$y, covariates = C1, n_iter = 20, burn_in = 10,
                   seed = 1, verbose = FALSE)
  expect_equal(fit1$covariates, c("(Intercept)", "x"))
  fit2 <- balexmed(sim$X, sim$m, sim$y, covariates = cbind(1, C1), n_iter = 20,
                   burn_in = 10, seed = 1, verbose = FALSE)
  expect_equal(fit2$k, 2)
  fit0 <- balexmed(sim$X, sim$m, sim$y, n_iter = 20, burn_in = 10, seed = 1,
                   sampler = "mh", verbose = FALSE)
  expect_equal(fit0$covariates, "(Intercept)")
})

test_that("balexmed validates its input", {
  set.seed(12)
  sim <- simulate_balexmed(n = 30, d = 6)
  X0 <- sim$X; X0[1, 1] <- 0
  expect_error(balexmed(X0, sim$m, sim$y, verbose = FALSE), "strictly positive")
  expect_error(balexmed(sim$X, sim$m[-1], sim$y, verbose = FALSE), "length 30")
  expect_error(balexmed(sim$X, sim$m, sim$y, eta = 0.5, verbose = FALSE), "eta")
  expect_error(balexmed(sim$X, sim$m, sim$y, n_iter = 10, burn_in = 10, verbose = FALSE),
               "No draws")
  expect_error(balexmed(sim$X, sim$m, sim$y, z_init = c(1, 1, 0, 0, 0, 0), verbose = FALSE),
               "at least one 1 and one -1")
  expect_error(balexmed(sim$X, sim$m, sim$y, covariates = data.frame(a = c(NA, 1:29)),
                        verbose = FALSE), "missing")
})
