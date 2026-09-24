## Collapsed posterior of the balance configuration z (internal).
##
## With the normal-inverse-gamma priors, the coefficients and error variances
## integrate out in closed form, so that
##   log f(z | y, m) = -1/2 log|Lambda_m| - (nu + n)/2 log(lambda_m + SS_m)
##                     -1/2 log|Lambda_y| - (nu + n)/2 log(lambda_y + SS_y)
##                     + log f(z) + constant,
## where Lambda_m = D_m'D_m + h I for D_m = (B_z, C) and Lambda_y = D_y'D_y + h I
## for D_y = (m, B_z, C). The covariate block C'C + h I does not depend on z, so
## it is factorized once and each evaluation reduces to Schur complements of
## size one (mediator model) and two (outcome model).

.bm_setup <- function(X, m, y, C, eta, h, nu, lambda_m, lambda_y) {
  G <- crossprod(C) + diag(h, ncol(C))
  R <- chol(G)
  um <- drop(backsolve(R, drop(crossprod(C, m)), transpose = TRUE))
  uy <- drop(backsolve(R, drop(crossprod(C, y)), transpose = TRUE))
  list(n = nrow(X), d = ncol(X), k = ncol(C), logX = log(X), C = C,
       m = m, y = y, R = R, logdet_G = 2 * sum(log(diag(R))),
       um = um, uy = uy,
       mtm = sum(m * m), yty = sum(y * y), mty = sum(m * y),
       um2 = sum(um * um), uy2 = sum(uy * uy), umuy = sum(um * uy),
       h = h, nu = nu, lambda_m = lambda_m, lambda_y = lambda_y,
       log_eta = log(eta), log_eta0 = log(1 - 2 * eta))
}

## Log posterior (up to a constant) for a balance vector B with ap numerator
## and am denominator parts.
.bm_logpost <- function(B, ap, am, s) {
  BtB <- sum(B * B)
  Btm <- sum(B * s$m)
  Bty <- sum(B * s$y)
  v <- drop(backsolve(s$R, drop(crossprod(s$C, B)), transpose = TRUE))
  vv <- sum(v * v)
  vum <- sum(v * s$um)
  vuy <- sum(v * s$uy)
  n <- s$n
  ## mediator model, design (B, C)
  schur_B <- BtB + s$h - vv
  ss_m <- s$mtm - s$um2 - (Btm - vum)^2 / schur_B
  ldet_m <- s$logdet_G + log(schur_B)
  ## outcome model, design (m, B, C); 2 x 2 Schur complement S
  S11 <- s$mtm + s$h - s$um2
  S12 <- Btm - vum
  S22 <- schur_B
  r1 <- s$mty - s$umuy
  r2 <- Bty - vuy
  det_S <- S11 * S22 - S12 * S12
  ss_y <- s$yty - s$uy2 - (S22 * r1 * r1 - 2 * S12 * r1 * r2 + S11 * r2 * r2) / det_S
  ldet_y <- s$logdet_G + log(det_S)
  -0.5 * ldet_m - 0.5 * (s$nu + n) * log(s$lambda_m + ss_m) -
    0.5 * ldet_y - 0.5 * (s$nu + n) * log(s$lambda_y + ss_y) +
    (ap + am) * s$log_eta + (s$d - ap - am) * s$log_eta0
}

## Draw (alpha, psi_m, sigma2_m) and (beta, gamma, psi_y, sigma2_y) given z,
## then relabel so that gamma > 0 (Algorithm 1 of the paper).
.bm_draw_coefficients <- function(z, s) {
  B <- .balance_from_log(s$logX, z)
  shape <- (s$nu + s$n) / 2

  Dm <- cbind(B, s$C)
  Rm <- chol(crossprod(Dm) + diag(s$h, ncol(Dm)))
  rhs_m <- drop(crossprod(Dm, s$m))
  mu_m <- drop(backsolve(Rm, backsolve(Rm, rhs_m, transpose = TRUE)))
  ss_m <- s$mtm - sum(rhs_m * mu_m)
  sigma2_m <- 1 / stats::rgamma(1, shape = shape, rate = (s$lambda_m + ss_m) / 2)
  coef_m <- mu_m + sqrt(sigma2_m) * drop(backsolve(Rm, stats::rnorm(ncol(Dm))))

  Dy <- cbind(s$m, B, s$C)
  Ry <- chol(crossprod(Dy) + diag(s$h, ncol(Dy)))
  rhs_y <- drop(crossprod(Dy, s$y))
  mu_y <- drop(backsolve(Ry, backsolve(Ry, rhs_y, transpose = TRUE)))
  ss_y <- s$yty - sum(rhs_y * mu_y)
  sigma2_y <- 1 / stats::rgamma(1, shape = shape, rate = (s$lambda_y + ss_y) / 2)
  coef_y <- mu_y + sqrt(sigma2_y) * drop(backsolve(Ry, stats::rnorm(ncol(Dy))))

  flip <- coef_y[2] < 0
  if (flip) {
    z <- -z
    coef_m[1] <- -coef_m[1]
    coef_y[2] <- -coef_y[2]
  }
  list(z = z, alpha = coef_m[1], psi_m = coef_m[-1],
       beta = coef_y[1], gamma = coef_y[2], psi_y = coef_y[-(1:2)],
       sigma2_m = sigma2_m, sigma2_y = sigma2_y, flip = flip)
}
