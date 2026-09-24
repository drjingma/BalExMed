test_that("balance matches its definition", {
  X <- matrix(c(0.1, 0.2, 0.3, 0.4,
                0.4, 0.3, 0.2, 0.1,
                0.25, 0.25, 0.25, 0.25), nrow = 3, byrow = TRUE)
  z <- c(1, 1, -1, 0)
  manual <- sqrt(2 * 1 / 3) * (rowMeans(log(X[, 1:2])) - log(X[, 3]))
  expect_equal(balance(X, z), manual)
})

test_that("balance is invariant to rescaling rows and antisymmetric in z", {
  set.seed(2)
  X <- matrix(rexp(40), 8, 5)
  z <- c(1, -1, 0, 1, -1)
  expect_equal(balance(X * runif(8, 1, 100), z), balance(X, z))
  expect_equal(balance(X, -z), -balance(X, z))
})

test_that("balance rejects invalid input", {
  X <- matrix(c(0.5, 0.5, 0, 1), 2)
  expect_error(balance(X, c(1, -1)), "strictly positive")
  X <- matrix(runif(6), 2)
  expect_error(balance(X, c(1, 1, 0)), "at least one 1 and one -1")
  expect_error(balance(X, c(1, -1)), "length 3")
  expect_error(balance(X, c(1, -1, 2)), "only -1, 0 and 1")
})

test_that("pca_start returns a valid configuration", {
  set.seed(3)
  sim <- simulate_balexmed(n = 60, d = 15)
  z <- pca_start(sim$X)
  expect_length(z, 15)
  expect_true(any(z == 1) && any(z == -1))
  expect_error(pca_start(sim$X, threshold = 0.99), "Lower the threshold")
})
