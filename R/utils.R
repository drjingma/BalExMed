## Input checks and small helpers (internal).

.as_composition <- function(X) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X))
    stop("`X` must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X))
    stop("`X` must not contain missing values.", call. = FALSE)
  if (any(!is.finite(X)) || any(X <= 0))
    stop("`X` must be strictly positive. Replace zeros before fitting, ",
         "for example with a pseudocount.", call. = FALSE)
  X
}

.check_z <- function(z, d, arg = "z") {
  if (!is.numeric(z) || length(z) != d)
    stop(sprintf("`%s` must be a numeric vector of length %d.", arg, d),
         call. = FALSE)
  if (anyNA(z) || !all(z %in% c(-1, 0, 1)))
    stop(sprintf("`%s` must contain only -1, 0 and 1.", arg), call. = FALSE)
  if (!any(z == 1) || !any(z == -1))
    stop(sprintf("`%s` must contain at least one 1 and one -1.", arg),
         call. = FALSE)
  as.numeric(z)
}

.as_response <- function(v, n, arg) {
  if (is.data.frame(v)) v <- as.matrix(v)
  if (is.matrix(v)) {
    if (ncol(v) != 1L)
      stop(sprintf("`%s` must be a vector or a one-column matrix.", arg),
           call. = FALSE)
    v <- v[, 1]
  }
  if (!is.numeric(v) || length(v) != n)
    stop(sprintf("`%s` must be a numeric vector of length %d, one value per row of `X`.",
                 arg, n), call. = FALSE)
  if (anyNA(v) || any(!is.finite(v)))
    stop(sprintf("`%s` must not contain missing or infinite values.", arg),
         call. = FALSE)
  as.numeric(v)
}

## Covariate design matrix. Always contains an intercept.
.covariate_matrix <- function(covariates, n) {
  if (is.null(covariates))
    return(matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)")))
  if (anyNA(covariates))
    stop("`covariates` must not contain missing values.", call. = FALSE)
  if (is.data.frame(covariates)) {
    if (nrow(covariates) != n)
      stop("`covariates` must have one row per row of `X`.", call. = FALSE)
    C <- stats::model.matrix(~ ., data = covariates)
    attr(C, "assign") <- NULL
    attr(C, "contrasts") <- NULL
  } else {
    C <- as.matrix(covariates)
    if (!is.numeric(C))
      stop("A covariate matrix must be numeric; pass a data frame to use factors.",
           call. = FALSE)
    if (nrow(C) != n)
      stop("`covariates` must have one row per row of `X`.", call. = FALSE)
    if (is.null(colnames(C))) colnames(C) <- paste0("covariate", seq_len(ncol(C)))
    constant <- apply(C, 2, function(v) all(v == v[1]))
    if (!any(constant)) C <- cbind(`(Intercept)` = 1, C)
  }
  if (any(!is.finite(C)))
    stop("`covariates` must not contain infinite values.", call. = FALSE)
  if (qr(C)$rank < ncol(C))
    stop("The covariate design matrix, including the intercept, is rank deficient.",
         call. = FALSE)
  C
}

.check_count <- function(x, arg, min = 1) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min || x != round(x))
    stop(sprintf("`%s` must be a whole number of at least %d.", arg, min),
         call. = FALSE)
  as.integer(x)
}

.check_positive <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0)
    stop(sprintf("`%s` must be a single positive number.", arg), call. = FALSE)
  x
}

.check_fit <- function(object) {
  if (!inherits(object, "balexmed"))
    stop("`object` must be a fit returned by balexmed().", call. = FALSE)
  invisible(object)
}

## Save and restore the global random number generator state.
.get_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  else NULL
}

.set_rng_state <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
}
