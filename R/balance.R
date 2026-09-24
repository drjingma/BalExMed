#' Compute a balance
#'
#' Computes, for each row of a composition, the balance
#' \deqn{B(z, x_i) = \sqrt{\frac{a_+ a_-}{a_+ + a_-}}
#'   \left( \frac{1}{a_+} \sum_{j: z_j = 1} \log x_{ij} -
#'          \frac{1}{a_-} \sum_{j: z_j = -1} \log x_{ij} \right),}
#' where \eqn{a_+} and \eqn{a_-} are the numbers of parts in the numerator and
#' the denominator.
#'
#' @param X An \eqn{n \times d} numeric matrix or data frame of strictly
#'   positive abundances; rows are samples and columns are parts. Rows need not
#'   sum to one because a balance does not change when a row is rescaled.
#' @param z A vector of length \eqn{d} with entries in \{-1, 0, 1\}. A value of
#'   1 places a part in the numerator, -1 in the denominator and 0 leaves it
#'   out. Both groups must be non-empty.
#' @return A numeric vector of length \eqn{n}.
#' @examples
#' X <- matrix(c(0.2, 0.3, 0.5,
#'               0.6, 0.1, 0.3), nrow = 2, byrow = TRUE)
#' balance(X, c(1, -1, 0))
#' sqrt(1 / 2) * log(X[, 1] / X[, 2])
#' @export
balance <- function(X, z) {
  X <- .as_composition(X)
  z <- .check_z(z, ncol(X))
  .balance_from_log(log(X), z)
}

.balance_from_log <- function(logX, z) {
  pos <- which(z == 1)
  neg <- which(z == -1)
  ap <- length(pos)
  am <- length(neg)
  sqrt(ap * am / (ap + am)) *
    (rowMeans(logX[, pos, drop = FALSE]) - rowMeans(logX[, neg, drop = FALSE]))
}

#' Principal-component starting value for the balance configuration
#'
#' Builds a starting value for \eqn{z} from the leading right singular vector
#' of the centered log-ratio transformed composition. Parts whose loading
#' exceeds `threshold` in absolute value enter the numerator or the
#' denominator according to the sign of the loading.
#'
#' @inheritParams balance
#' @param threshold Absolute loading above which a part is included. The
#'   loading vector has unit length, so the default of 0.1 suits compositions
#'   with tens to a few hundred parts.
#' @return A vector with entries in \{-1, 0, 1\}.
#' @examples
#' set.seed(1)
#' sim <- simulate_balexmed(n = 50, d = 20)
#' pca_start(sim$X)
#' @export
pca_start <- function(X, threshold = 0.1) {
  X <- .as_composition(X)
  threshold <- .check_positive(threshold, "threshold")
  lx <- log(X)
  clr <- lx - rowMeans(lx)
  v <- svd(clr, nu = 0, nv = 1)$v[, 1]
  z <- sign(v) * (abs(v) > threshold)
  if (!any(z == 1) || !any(z == -1))
    stop("The leading loading vector has no parts of one sign above `threshold`. ",
         "Lower the threshold or use a random start.", call. = FALSE)
  z
}

.random_start <- function(d, prob_active = 0.1) {
  z <- numeric(d)
  k <- max(2L, ceiling(prob_active * d))
  active <- sample.int(d, size = k)
  z[active] <- sample(c(-1, 1), size = k, replace = TRUE)
  if (!any(z == 1)) z[active[sample.int(k, 1L)]] <- 1
  if (!any(z == -1)) z[active[sample.int(k, 1L)]] <- -1
  z
}
