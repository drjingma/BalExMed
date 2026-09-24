#' Fit the balance mediation model for a compositional exposure
#'
#' Estimates the balance mediation model
#' \deqn{m_i = \alpha B(z, x_i) + \psi_m^\top c_i + e_i, \qquad
#'       y_i = \gamma B(z, x_i) + \beta m_i + \psi_y^\top c_i + \varepsilon_i,}
#' where \eqn{B(z, x_i)} is the balance defined in [balance()], \eqn{e_i \sim
#' N(0, \sigma_e^2)} and \eqn{\varepsilon_i \sim N(0, \sigma_\varepsilon^2)}.
#' The configuration \eqn{z} is latent and is sampled jointly with the path
#' coefficients.
#'
#' @section Priors:
#' The coefficients have conjugate normal priors whose variances are the error
#' variance divided by `h`, and each error variance has an inverse-gamma prior
#' with shape `nu / 2` and rate `lambda_m / 2` or `lambda_y / 2`. Each part
#' independently enters the numerator with probability `eta`, the denominator
#' with probability `eta`, and neither with probability `1 - 2 * eta`, so the
#' expected number of parts on each side is `eta * d`. The value `eta = 1/3` is
#' non-informative; smaller values give sparser balances.
#'
#' @section Sampling:
#' Coefficients and variances are integrated out, so \eqn{z} is updated from its
#' collapsed posterior. The Gibbs sampler updates every coordinate of \eqn{z}
#' once per iteration. The Metropolis-Hastings sampler changes at most one
#' coordinate per iteration; it is cheaper per iteration but mixes more slowly.
#' After \eqn{z} is updated, the coefficients and variances are drawn from their
#' conditional posteriors. Exchanging the numerator and denominator changes the
#' sign of the balance, so each draw is relabeled to satisfy \eqn{\gamma > 0}:
#' the numerator holds the parts whose higher relative abundance is associated
#' with a higher outcome through the direct path.
#'
#' Zeros must be replaced before fitting because the balance uses logarithms.
#' Standardizing continuous covariates is recommended but not required.
#'
#' @param X An \eqn{n \times d} numeric matrix or data frame of strictly
#'   positive abundances, with samples in rows and parts (for example, taxa) in
#'   columns. Rows are rescaled to sum to one. Column names are used as part
#'   names.
#' @param m Numeric mediator, one value per row of `X`.
#' @param y Numeric outcome, one value per row of `X`.
#' @param covariates Optional covariates adjusted for in both regressions: a
#'   data frame, expanded with [stats::model.matrix()] so factors are allowed,
#'   or a numeric matrix. An intercept is always included.
#' @param eta Prior probability that a part enters the numerator, and also the
#'   prior probability that it enters the denominator. Must lie in (0, 1/2).
#' @param n_iter Total number of iterations per chain, including burn-in.
#' @param burn_in Number of initial iterations discarded from each chain.
#' @param thin Keep every `thin`-th iteration after burn-in.
#' @param n_chains Number of independent chains.
#' @param sampler Either `"gibbs"` (default) or `"mh"` for Metropolis-Hastings.
#' @param z_init Starting configuration: `NULL` for random starts in every
#'   chain; `"pca"` to start the first chain at [pca_start()] and the others at
#'   random; a vector used by every chain; or a list with one vector (or
#'   `NULL`) per chain.
#' @param seed Optional integer. Chain `i` is seeded with `seed + i - 1`, so
#'   results do not depend on `cores`. When `NULL`, chain seeds are drawn from
#'   the current random number stream. The global random number state is left
#'   as it was after those seeds were drawn.
#' @param cores Number of cores used to run chains in parallel with
#'   [parallel::mclapply()]. Parallel execution is not available on Windows,
#'   where chains run sequentially.
#' @param h Prior precision factor of every regression coefficient.
#' @param nu Prior degrees of freedom of both error variances; defaults to `n`.
#' @param lambda_m,lambda_y Prior scale of the mediator and outcome error
#'   variances; default to the sample variances of `m` and `y`.
#' @param verbose Print progress messages.
#'
#' @return An object of class `"balexmed"`: a list with elements `chains` (one
#'   list of retained draws per chain), `taxa`, `covariates` (design column
#'   names), `n`, `d`, `k` and `settings`. Use [summary.balexmed()],
#'   [taxa_table()], [effect_draws()], [selected_z()] and [plot.balexmed()] to
#'   inspect it.
#'
#' @seealso [balance()], [simulate_balexmed()]
#' @examples
#' set.seed(1)
#' sim <- simulate_balexmed(n = 100, d = 20)
#' fit <- balexmed(sim$X, sim$m, sim$y, covariates = sim$covariates,
#'                 n_iter = 400, burn_in = 200, seed = 1, verbose = FALSE)
#' summary(fit)
#' head(taxa_table(fit))
#' @export
balexmed <- function(X, m, y, covariates = NULL, eta = 1/5,
                     n_iter = 20000, burn_in = floor(n_iter / 2), thin = 1,
                     n_chains = 1, sampler = c("gibbs", "mh"), z_init = NULL,
                     seed = NULL, cores = 1, h = 1e-6, nu = NULL,
                     lambda_m = NULL, lambda_y = NULL, verbose = TRUE) {
  cl <- match.call()
  sampler <- match.arg(sampler)
  X <- .as_composition(X)
  n <- nrow(X)
  d <- ncol(X)
  if (d < 2L) stop("`X` must have at least two columns.", call. = FALSE)
  taxa <- colnames(X)
  if (is.null(taxa)) taxa <- paste0("part", seq_len(d))
  X <- X / rowSums(X)
  m <- .as_response(m, n, "m")
  y <- .as_response(y, n, "y")
  C <- .covariate_matrix(covariates, n)
  if (n <= ncol(C) + 2L)
    stop("The number of samples must exceed the number of covariate columns plus two.",
         call. = FALSE)

  if (!is.numeric(eta) || length(eta) != 1L || is.na(eta) || eta <= 0 || eta >= 0.5)
    stop("`eta` must be a single number in (0, 1/2).", call. = FALSE)
  n_iter <- .check_count(n_iter, "n_iter")
  burn_in <- .check_count(burn_in, "burn_in", min = 0)
  thin <- .check_count(thin, "thin")
  n_chains <- .check_count(n_chains, "n_chains")
  cores <- .check_count(cores, "cores")
  if ((n_iter - burn_in) %/% thin < 1L)
    stop("No draws would be kept: `n_iter` must exceed `burn_in` by at least `thin`.",
         call. = FALSE)
  h <- .check_positive(h, "h")
  nu <- if (is.null(nu)) n else .check_positive(nu, "nu")
  lambda_m <- if (is.null(lambda_m)) stats::var(m) else .check_positive(lambda_m, "lambda_m")
  lambda_y <- if (is.null(lambda_y)) stats::var(y) else .check_positive(lambda_y, "lambda_y")
  if (lambda_m <= 0 || lambda_y <= 0)
    stop("`m` and `y` must vary, or `lambda_m` and `lambda_y` must be supplied.",
         call. = FALSE)

  s <- .bm_setup(X, m, y, C, eta, h, nu, lambda_m, lambda_y)
  inits <- .resolve_inits(z_init, X, n_chains)
  seeds <- .chain_seeds(seed, n_chains)
  rng_after_seeds <- .get_rng_state()
  on.exit(.set_rng_state(rng_after_seeds), add = TRUE)

  run_chain <- function(i) {
    set.seed(seeds[i])
    z0 <- if (is.null(inits[[i]])) .random_start(d) else inits[[i]]
    t0 <- proc.time()[["elapsed"]]
    res <- .bm_chain(s, n_iter, burn_in, thin, z0, sampler)
    res$seconds <- proc.time()[["elapsed"]] - t0
    res$seed <- seeds[i]
    res$z_init <- stats::setNames(z0, taxa)
    colnames(res$z) <- taxa
    colnames(res$psi_m) <- colnames(C)
    colnames(res$psi_y) <- colnames(C)
    res
  }

  if (verbose)
    message(sprintf("Fitting %d chain(s) with the %s sampler: n = %d, d = %d, %d covariate column(s), %d iterations.",
                    n_chains, if (sampler == "gibbs") "Gibbs" else "Metropolis-Hastings",
                    n, d, ncol(C), n_iter))
  use_fork <- cores > 1L && n_chains > 1L && .Platform$OS.type != "windows"
  chains <- if (use_fork) {
    parallel::mclapply(seq_len(n_chains), run_chain,
                       mc.cores = min(cores, n_chains), mc.set.seed = FALSE)
  } else {
    lapply(seq_len(n_chains), run_chain)
  }
  failed <- vapply(chains, function(ch) inherits(ch, "try-error") || is.null(ch$gamma),
                   logical(1))
  if (any(failed))
    stop("Chain(s) ", paste(which(failed), collapse = ", "), " failed: ",
         paste(unique(vapply(chains[failed], as.character, character(1))), collapse = "; "),
         call. = FALSE)
  if (verbose)
    message(sprintf("Done in %.1f seconds of chain time.",
                    sum(vapply(chains, `[[`, numeric(1), "seconds"))))

  structure(list(
    chains = chains,
    call = cl,
    taxa = taxa,
    covariates = colnames(C),
    n = n, d = d, k = ncol(C),
    settings = list(sampler = sampler, eta = eta, h = h, nu = nu,
                    lambda_m = lambda_m, lambda_y = lambda_y,
                    n_iter = n_iter, burn_in = burn_in, thin = thin,
                    n_chains = n_chains)
  ), class = "balexmed")
}

.resolve_inits <- function(z_init, X, n_chains) {
  d <- ncol(X)
  out <- vector("list", n_chains)
  if (is.null(z_init)) return(out)
  if (is.character(z_init)) {
    if (!identical(z_init, "pca"))
      stop("`z_init` must be NULL, \"pca\", a vector, or a list of vectors.", call. = FALSE)
    out[[1]] <- pca_start(X)
    return(out)
  }
  if (is.list(z_init)) {
    if (length(z_init) != n_chains)
      stop("A list `z_init` must have one element per chain.", call. = FALSE)
    for (i in seq_len(n_chains))
      if (!is.null(z_init[[i]])) out[[i]] <- .check_z(z_init[[i]], d, "z_init")
    return(out)
  }
  v <- .check_z(z_init, d, "z_init")
  for (i in seq_len(n_chains)) out[[i]] <- v
  out
}

.chain_seeds <- function(seed, n_chains) {
  if (is.null(seed)) return(sample.int(.Machine$integer.max, n_chains))
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed != round(seed))
    stop("`seed` must be a single whole number or NULL.", call. = FALSE)
  as.integer((seed + seq_len(n_chains) - 1) %% .Machine$integer.max)
}
