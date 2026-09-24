## Updates of z and the full MCMC chain (internal).

## One Gibbs sweep: visit the d coordinates in random order and draw each z_j
## from its three-category full conditional, restricted to configurations with
## non-empty numerator and denominator.
.bm_gibbs_sweep <- function(z, s) {
  logX <- s$logX
  pos <- z == 1
  neg <- z == -1
  Sp <- rowSums(logX[, pos, drop = FALSE])
  Sm <- rowSums(logX[, neg, drop = FALSE])
  ap <- sum(pos)
  am <- sum(neg)
  n_moved <- 0L
  states <- c(-1, 0, 1)
  for (j in sample.int(s$d, size = s$d)) {
    zj <- z[j]
    xj <- logX[, j]
    Sp0 <- if (zj == 1) Sp - xj else Sp
    ap0 <- ap - (zj == 1)
    Sm0 <- if (zj == -1) Sm - xj else Sm
    am0 <- am - (zj == -1)
    lp <- c(-Inf, -Inf, -Inf)
    if (ap0 > 0) {                       # z_j = -1
      a <- am0 + 1
      B <- sqrt(ap0 * a / (ap0 + a)) * (Sp0 / ap0 - (Sm0 + xj) / a)
      lp[1] <- .bm_logpost(B, ap0, a, s)
    }
    if (ap0 > 0 && am0 > 0) {            # z_j = 0
      B <- sqrt(ap0 * am0 / (ap0 + am0)) * (Sp0 / ap0 - Sm0 / am0)
      lp[2] <- .bm_logpost(B, ap0, am0, s)
    }
    if (am0 > 0) {                       # z_j = 1
      a <- ap0 + 1
      B <- sqrt(a * am0 / (a + am0)) * ((Sp0 + xj) / a - Sm0 / am0)
      lp[3] <- .bm_logpost(B, a, am0, s)
    }
    ok <- is.finite(lp)
    if (!any(ok)) next
    prob <- exp(lp - max(lp[ok]))
    prob[!ok] <- 0
    new <- states[sample.int(3L, 1L, prob = prob / sum(prob))]
    if (new != zj) {
      n_moved <- n_moved + 1L
      z[j] <- new
      Sp <- if (new == 1) Sp0 + xj else Sp0
      ap <- ap0 + (new == 1)
      Sm <- if (new == -1) Sm0 + xj else Sm0
      am <- am0 + (new == -1)
    }
  }
  list(z = z, n_moved = n_moved)
}

## One Metropolis-Hastings step: pick one of the value pairs (0, 1), (0, -1),
## (1, -1) with equal probability, pick a coordinate uniformly among those
## taking a value in the pair, and propose the other value of the pair. The
## number of eligible coordinates is the same before and after the move, so
## the proposal is symmetric.
.bm_mh_step <- function(z, lp, s) {
  pairs <- list(c(0, 1), c(0, -1), c(1, -1))
  pr <- pairs[[sample.int(3L, 1L)]]
  elig <- which(z == pr[1] | z == pr[2])
  if (length(elig) == 0L) return(list(z = z, lp = lp, accepted = FALSE))
  j <- elig[sample.int(length(elig), 1L)]
  z_new <- z
  z_new[j] <- if (z[j] == pr[1]) pr[2] else pr[1]
  ap <- sum(z_new == 1)
  am <- sum(z_new == -1)
  if (ap == 0 || am == 0) return(list(z = z, lp = lp, accepted = FALSE))
  lp_new <- .bm_logpost(.balance_from_log(s$logX, z_new), ap, am, s)
  if (log(stats::runif(1)) < lp_new - lp)
    return(list(z = z_new, lp = lp_new, accepted = TRUE))
  list(z = z, lp = lp, accepted = FALSE)
}

.bm_chain <- function(s, n_iter, burn_in, thin, z_init, sampler) {
  n_keep <- (n_iter - burn_in) %/% thin
  d <- s$d
  k <- s$k
  z_draws <- matrix(0L, n_keep, d)
  alpha <- beta <- gamma <- sigma2_m <- sigma2_y <- numeric(n_keep)
  psi_m <- psi_y <- matrix(0, n_keep, k)
  moves <- 0
  flips <- 0
  z <- z_init
  lp <- if (sampler == "mh")
    .bm_logpost(.balance_from_log(s$logX, z), sum(z == 1), sum(z == -1), s)
  else NA_real_
  idx <- 0L
  for (it in seq_len(n_iter)) {
    if (sampler == "gibbs") {
      sw <- .bm_gibbs_sweep(z, s)
      z <- sw$z
      moves <- moves + (sw$n_moved > 0L)
    } else {
      st <- .bm_mh_step(z, lp, s)
      z <- st$z
      lp <- st$lp
      moves <- moves + st$accepted
    }
    ## the posterior of z is symmetric under z -> -z, so lp is unchanged by a flip
    cf <- .bm_draw_coefficients(z, s)
    z <- cf$z
    flips <- flips + cf$flip
    if (it > burn_in && (it - burn_in) %% thin == 0L && idx < n_keep) {
      idx <- idx + 1L
      z_draws[idx, ] <- as.integer(z)
      alpha[idx] <- cf$alpha
      beta[idx] <- cf$beta
      gamma[idx] <- cf$gamma
      sigma2_m[idx] <- cf$sigma2_m
      sigma2_y[idx] <- cf$sigma2_y
      psi_m[idx, ] <- cf$psi_m
      psi_y[idx, ] <- cf$psi_y
    }
  }
  list(z = z_draws, alpha = alpha, beta = beta, gamma = gamma,
       sigma2_m = sigma2_m, sigma2_y = sigma2_y, psi_m = psi_m, psi_y = psi_y,
       move_rate = moves / n_iter, flip_rate = flips / n_iter)
}
