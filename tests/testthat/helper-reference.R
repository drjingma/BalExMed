## Reference implementations used only by the tests.

## Collapsed log posterior of z computed directly from the (k+1)- and
## (k+2)-dimensional precision matrices, without the Schur-complement shortcuts.
naive_logpost <- function(z, X, m, y, C, eta, h, nu, lambda_m, lambda_y) {
  n <- nrow(X)
  k <- ncol(C)
  B <- balance(X, z)
  Dm <- cbind(B, C)
  Lm <- crossprod(Dm) + diag(h, k + 1)
  mu_m <- solve(Lm, crossprod(Dm, m))
  ss_m <- sum(m^2) - drop(t(mu_m) %*% Lm %*% mu_m)
  Dy <- cbind(m, B, C)
  Ly <- crossprod(Dy) + diag(h, k + 2)
  mu_y <- solve(Ly, crossprod(Dy, y))
  ss_y <- sum(y^2) - drop(t(mu_y) %*% Ly %*% mu_y)
  prior <- sum(z != 0) * log(eta) + sum(z == 0) * log(1 - 2 * eta)
  -0.5 * determinant(Lm)$modulus[1] - (nu + n) / 2 * log(lambda_m + ss_m) -
    0.5 * determinant(Ly)$modulus[1] - (nu + n) / 2 * log(lambda_y + ss_y) +
    prior
}

## Exact log marginal density of (m, y) given z, computed with n x n matrices:
## v | z follows a multivariate t with scale matrix I + D H^{-1} D'.
exact_logpost <- function(z, X, m, y, C, eta, h, nu, lambda_m, lambda_y) {
  n <- nrow(X)
  B <- balance(X, z)
  block <- function(v, D, lambda) {
    S <- diag(n) + tcrossprod(D) / h
    -0.5 * determinant(S)$modulus[1] -
      (nu + n) / 2 * log(lambda + drop(t(v) %*% solve(S, v)))
  }
  prior <- sum(z != 0) * log(eta) + sum(z == 0) * log(1 - 2 * eta)
  block(m, cbind(B, C), lambda_m) + block(y, cbind(m, B, C), lambda_y) + prior
}

## Posterior probability that gamma > 0 given z (marginal t distribution).
prob_gamma_positive <- function(z, X, m, y, C, h, nu, lambda_y) {
  n <- nrow(X)
  B <- balance(X, z)
  Dy <- cbind(m, B, C)
  Ly <- crossprod(Dy) + diag(h, ncol(Dy))
  Li <- solve(Ly)
  mu <- drop(Li %*% crossprod(Dy, y))
  ss <- sum(y^2) - sum(mu * drop(crossprod(Dy, y)))
  scale <- sqrt((lambda_y + ss) / (nu + n) * Li[2, 2])
  stats::pt(mu[2] / scale, df = nu + n)
}

all_valid_z <- function(d) {
  g <- as.matrix(expand.grid(rep(list(c(-1, 0, 1)), d)))
  g[apply(g, 1, function(z) any(z == 1) && any(z == -1)), , drop = FALSE]
}
