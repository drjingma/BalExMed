#' Simulate data from the balance mediation model
#'
#' Generates data following the simulation design of the accompanying paper.
#' Log-abundances are multivariate normal with covariance
#' \eqn{\Sigma_{jj'} = v \rho^{|j - j'|}}; they are exponentiated and rounded
#' to counts, zeros are replaced by 0.5, and each row is closed to sum to one.
#' The first `n_plus` parts form the numerator and the next `n_minus` parts the
#' denominator. One standard normal covariate enters both regressions with
#' coefficient 1, and both error terms are standard normal.
#'
#' @param n Number of samples.
#' @param d Number of parts.
#' @param alpha,beta,gamma Path coefficients: balance to mediator, mediator to
#'   outcome, and balance to outcome.
#' @param n_plus,n_minus Numbers of parts in the numerator and denominator.
#' @param rho Correlation parameter of the log-abundances.
#' @param var_log Variance \eqn{v} of each log-abundance.
#' @return A list with `X` (compositions), `m`, `y`, `covariates` (a data
#'   frame), the true configuration `z`, the true `balance`, and `effects`
#'   (the true direct and indirect effects).
#' @examples
#' set.seed(1)
#' sim <- simulate_balexmed(n = 100, d = 30)
#' str(sim, max.level = 1)
#' @export
simulate_balexmed <- function(n = 100, d = 50, alpha = 1/3, beta = 2, gamma = 1/3,
                              n_plus = 3, n_minus = 3, rho = 0.2, var_log = 9) {
  n <- .check_count(n, "n", min = 3)
  d <- .check_count(d, "d", min = 2)
  n_plus <- .check_count(n_plus, "n_plus")
  n_minus <- .check_count(n_minus, "n_minus")
  if (n_plus + n_minus > d)
    stop("`n_plus + n_minus` must not exceed `d`.", call. = FALSE)
  S <- var_log * rho^abs(outer(seq_len(d), seq_len(d), "-"))
  log_abund <- matrix(stats::rnorm(n * d), n, d) %*% chol(S)
  counts <- round(exp(log_abund))
  counts[counts == 0] <- 0.5
  X <- counts / rowSums(counts)
  colnames(X) <- paste0("taxon", seq_len(d))
  z <- c(rep(1, n_plus), rep(-1, n_minus), rep(0, d - n_plus - n_minus))
  B <- balance(X, z)
  covariate <- stats::rnorm(n)
  m <- alpha * B + covariate + stats::rnorm(n)
  y <- gamma * B + beta * m + covariate + stats::rnorm(n)
  list(X = X, m = m, y = y, covariates = data.frame(c = covariate),
       z = stats::setNames(z, colnames(X)), balance = B,
       effects = c(direct = gamma, indirect = alpha * beta))
}
