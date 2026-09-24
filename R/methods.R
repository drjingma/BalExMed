#' Posterior draws of the path coefficients and effects
#'
#' @param object A fit returned by [balexmed()].
#' @return A data frame with one row per retained draw and columns `chain`,
#'   `draw`, `alpha`, `beta`, `gamma` (direct effect), `indirect`
#'   (\eqn{\alpha\beta}) and `total` (\eqn{\gamma + \alpha\beta}).
#' @export
effect_draws <- function(object) {
  .check_fit(object)
  do.call(rbind, lapply(seq_along(object$chains), function(i) {
    ch <- object$chains[[i]]
    data.frame(chain = i, draw = seq_along(ch$alpha),
               alpha = ch$alpha, beta = ch$beta, gamma = ch$gamma,
               indirect = ch$alpha * ch$beta,
               total = ch$gamma + ch$alpha * ch$beta)
  }))
}

#' Posterior inclusion probabilities of the parts
#'
#' @param object A fit returned by [balexmed()].
#' @param threshold Posterior inclusion probability above which a part is
#'   flagged as selected.
#' @return A data frame with one row per part, ordered by decreasing posterior
#'   inclusion probability, with columns `taxon`; `P_plus` and `P_minus`, the
#'   posterior probabilities of the numerator and the denominator; `PIP`, their
#'   sum; `signed`, `P_plus - P_minus`; `side`; `selected`; and, for more than
#'   one chain, `PIP_sd_chains`, the standard deviation of the PIP across chains.
#' @export
taxa_table <- function(object, threshold = 0.5) {
  .check_fit(object)
  d <- object$d
  p_plus <- matrix(vapply(object$chains, function(ch) colMeans(ch$z == 1L),
                          numeric(d)), nrow = d)
  p_minus <- matrix(vapply(object$chains, function(ch) colMeans(ch$z == -1L),
                           numeric(d)), nrow = d)
  pip <- p_plus + p_minus
  out <- data.frame(taxon = object$taxa,
                    P_plus = rowMeans(p_plus),
                    P_minus = rowMeans(p_minus),
                    PIP = rowMeans(pip),
                    stringsAsFactors = FALSE)
  out$signed <- out$P_plus - out$P_minus
  out$side <- ifelse(out$signed > 0, "numerator",
                     ifelse(out$signed < 0, "denominator", "none"))
  out$selected <- out$PIP > threshold
  if (ncol(pip) > 1L) out$PIP_sd_chains <- apply(pip, 1, stats::sd)
  out <- out[order(-out$PIP, -abs(out$signed)), ]
  rownames(out) <- NULL
  out
}

#' Selected balance configuration
#'
#' Returns the configuration obtained by thresholding the posterior inclusion
#' probabilities: a part with PIP above `threshold` is placed on the side with
#' the larger posterior probability, and all other parts are set to 0. The
#' result can be passed to [balance()].
#'
#' @inheritParams taxa_table
#' @return A named vector with entries in \{-1, 0, 1\}, in the column order of
#'   `X`.
#' @export
selected_z <- function(object, threshold = 0.5) {
  tt <- taxa_table(object, threshold)
  tt <- tt[match(object$taxa, tt$taxon), ]
  z <- ifelse(tt$selected, sign(tt$signed), 0)
  if (!any(z == 1) || !any(z == -1))
    warning("The selected configuration does not have both a numerator and a ",
            "denominator, so its balance is undefined.", call. = FALSE)
  stats::setNames(z, object$taxa)
}

#' Convert draws of the effects to a coda object
#'
#' @param object A fit returned by [balexmed()].
#' @param pars Quantities to include; any of `"alpha"`, `"beta"`, `"gamma"`,
#'   `"indirect"` and `"total"`.
#' @return A [coda::mcmc.list()] with one element per chain.
#' @export
as_mcmc_list <- function(object, pars = c("alpha", "beta", "gamma", "indirect", "total")) {
  .check_fit(object)
  pars <- match.arg(pars, several.ok = TRUE)
  dr <- effect_draws(object)
  st <- object$settings
  coda::as.mcmc.list(lapply(split(dr, dr$chain), function(dd)
    coda::mcmc(as.matrix(dd[, pars, drop = FALSE]),
               start = st$burn_in + st$thin, thin = st$thin)))
}

## split-Rhat: each chain is split in half and the halves are compared
.split_rhat <- function(chains) {
  n <- min(lengths(chains))
  half <- floor(n / 2)
  if (half < 2L) return(NA_real_)
  subs <- unlist(lapply(chains, function(x) list(x[seq_len(half)], x[half + seq_len(half)])),
                 recursive = FALSE)
  means <- vapply(subs, mean, numeric(1))
  vars <- vapply(subs, stats::var, numeric(1))
  W <- mean(vars)
  if (!is.finite(W) || W <= 0) return(NA_real_)
  B <- half * stats::var(means)
  sqrt(((half - 1) / half * W + B / half) / W)
}

.ess <- function(chains) {
  n <- min(lengths(chains))
  if (n < 10L) return(NA_real_)
  tryCatch(
    sum(vapply(chains, function(x) unname(coda::effectiveSize(coda::mcmc(x[seq_len(n)]))),
               numeric(1))),
    error = function(e) NA_real_)
}

#' Summarize a balance mediation fit
#'
#' @param object A fit returned by [balexmed()].
#' @param level Probability of the equal-tailed credible intervals.
#' @param threshold Posterior inclusion probability above which a part is
#'   counted as selected.
#' @param x An object of class `"summary.balexmed"`.
#' @param digits Number of significant digits to print.
#' @param ... Unused.
#' @return An object of class `"summary.balexmed"`. Its element `effects` is a
#'   matrix with the posterior mean, standard deviation, credible interval,
#'   posterior probability of a positive value, split-Rhat and effective sample
#'   size of \eqn{\alpha}, \eqn{\beta}, the direct effect \eqn{\gamma}, the
#'   indirect effect \eqn{\alpha\beta} and the total effect.
#' @export
summary.balexmed <- function(object, level = 0.95, threshold = 0.5, ...) {
  .check_fit(object)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be a number in (0, 1).", call. = FALSE)
  dr <- effect_draws(object)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  pars <- c("alpha", "beta", "gamma", "indirect", "total")
  effects <- t(vapply(pars, function(p) {
    v <- dr[[p]]
    by_chain <- split(v, dr$chain)
    c(mean(v), stats::sd(v), stats::quantile(v, probs, names = FALSE),
      mean(v > 0), .split_rhat(by_chain), .ess(by_chain))
  }, numeric(7)))
  colnames(effects) <- c("mean", "sd", sprintf("%g%%", 100 * probs), "P(>0)", "Rhat", "ESS")
  rownames(effects) <- c("alpha (balance -> mediator)",
                         "beta (mediator -> outcome)",
                         "gamma (direct effect)",
                         "alpha*beta (indirect effect)",
                         "total effect")
  tt <- taxa_table(object, threshold)
  pip_cor <- NA_real_
  if (length(object$chains) > 1L) {
    pips <- vapply(object$chains, function(ch) colMeans(ch$z != 0L), numeric(object$d))
    cc <- suppressWarnings(stats::cor(pips))
    pip_cor <- min(cc[upper.tri(cc)], na.rm = TRUE)
  }
  structure(list(
    effects = effects, n = object$n, d = object$d, k = object$k,
    settings = object$settings, level = level, threshold = threshold,
    n_numerator = sum(tt$selected & tt$side == "numerator"),
    n_denominator = sum(tt$selected & tt$side == "denominator"),
    move_rate = vapply(object$chains, `[[`, numeric(1), "move_rate"),
    min_pip_cor = pip_cor
  ), class = "summary.balexmed")
}

#' @rdname summary.balexmed
#' @export
print.summary.balexmed <- function(x, digits = 3, ...) {
  st <- x$settings
  cat("Balance mediation model for a compositional exposure\n")
  cat(sprintf("n = %d samples, d = %d parts, %d covariate column(s) including the intercept\n",
              x$n, x$d, x$k))
  cat(sprintf("%s sampler: %d chain(s), %d iterations, burn-in %d, thin %d; eta = %g\n",
              if (st$sampler == "gibbs") "Gibbs" else "Metropolis-Hastings",
              st$n_chains, st$n_iter, st$burn_in, st$thin, st$eta))
  cat(sprintf("%s: %s\n",
              if (st$sampler == "gibbs") "Move rate" else "Acceptance rate",
              paste(formatC(x$move_rate, digits = 3, format = "f"), collapse = ", ")))
  cat(sprintf("\nPosterior summaries (%g%% credible intervals):\n", 100 * x$level))
  print(signif(x$effects, digits))
  cat(sprintf("\nParts with PIP > %g: %d in the numerator, %d in the denominator\n",
              x$threshold, x$n_numerator, x$n_denominator))
  if (!is.na(x$min_pip_cor))
    cat(sprintf("Minimum pairwise correlation of PIPs across chains: %.3f\n", x$min_pip_cor))
  invisible(x)
}

#' @export
print.balexmed <- function(x, ...) {
  st <- x$settings
  cat("Balance mediation fit (class \"balexmed\")\n")
  cat(sprintf("n = %d, d = %d, %s sampler, %d chain(s) of %d iterations\n",
              x$n, x$d, st$sampler, st$n_chains, st$n_iter))
  dr <- effect_draws(x)
  cat(sprintf("Posterior means: direct effect %.4g, indirect effect %.4g\n",
              mean(dr$gamma), mean(dr$indirect)))
  cat("Use summary() for credible intervals and taxa_table() for inclusion probabilities.\n")
  invisible(x)
}

#' Plot a balance mediation fit
#'
#' @param x A fit returned by [balexmed()].
#' @param type `"trace"` draws trace plots and running means of the direct and
#'   indirect effects, one line per chain. `"inclusion"` draws the signed
#'   inclusion probabilities of the parts with the largest PIP.
#' @param top Number of parts shown when `type = "inclusion"`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
plot.balexmed <- function(x, type = c("trace", "inclusion"), top = 30, ...) {
  type <- match.arg(type)
  if (type == "trace") {
    de <- vapply(x$chains, `[[`, numeric(length(x$chains[[1]]$gamma)), "gamma")
    ie <- vapply(x$chains, function(ch) ch$alpha * ch$beta, numeric(length(x$chains[[1]]$gamma)))
    de <- matrix(de, ncol = length(x$chains))
    ie <- matrix(ie, ncol = length(x$chains))
    running <- function(M) apply(M, 2, function(v) cumsum(v) / seq_along(v))
    cols <- seq_len(ncol(de)) + 1L
    op <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
    on.exit(graphics::par(op))
    graphics::matplot(de, type = "l", lty = 1, col = cols, xlab = "retained draw",
                      ylab = "gamma", main = "Direct effect")
    graphics::matplot(ie, type = "l", lty = 1, col = cols, xlab = "retained draw",
                      ylab = "alpha * beta", main = "Indirect effect")
    graphics::abline(h = 0, lty = 3)
    graphics::matplot(running(de), type = "l", lty = 1, col = cols, xlab = "retained draw",
                      ylab = "running mean", main = "Direct effect, running mean")
    graphics::matplot(running(ie), type = "l", lty = 1, col = cols, xlab = "retained draw",
                      ylab = "running mean", main = "Indirect effect, running mean")
  } else {
    tt <- utils::head(taxa_table(x), top)
    tt <- tt[rev(seq_len(nrow(tt))), ]
    op <- graphics::par(mar = c(4, max(4, 0.45 * max(nchar(tt$taxon))), 2.5, 1))
    on.exit(graphics::par(op))
    graphics::barplot(tt$signed, names.arg = tt$taxon, horiz = TRUE, las = 1,
                      xlim = c(-1, 1), col = ifelse(tt$signed > 0, "grey30", "grey75"),
                      xlab = "P(numerator) - P(denominator)",
                      main = "Signed inclusion probability")
    graphics::abline(v = 0)
  }
  invisible(x)
}
